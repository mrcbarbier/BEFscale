from utilities import *
import scipy.integrate as scint

#===========================================================================================================

# MODEL

#===========================================================================================================

class LandscapeModel():
    """Main Model object: given parameters, will generate matrices, run dynamics, load and save."""

    #Default parameters
    dft_prm=eval('\n'.join(line for line in open('default_parameters.dat')) ,{},{})

    methods=('scipy','Euler') #Available algorithms

    def __init__(self,**kwargs):
        """When the model object is created, store passed information."""
        self.data=deepcopy(kwargs.pop('data',{}))
        self.results=deepcopy(kwargs.pop('results',{}))
        self.prm=deepcopy(kwargs.pop('parameters',LandscapeModel.dft_prm))
        self.set_params(**kwargs)
        self.PARALLEL_LOCK=0


    def save(self,path,suffix=''):
        """Save model parameters in prm.dat, dynamical variables in results.npz and other matrices in data.npz"""
        fpath=Path(path)
        if not self.PARALLEL_LOCK:
            # Prevent overwriting data if parallel processing
            fpath.mkdir()
            dumps(open(fpath+'prm.dat','w'),self.prm)
            np.savez(fpath+'data',**self.data)
        np.savez(fpath+'results{}'.format(suffix),**self.results)


    @classmethod
    def load(cls,path):
        """Load variables and other parameters from save files in path."""
        fpath=Path(path)
        prm=loads(open(fpath+'prm.dat','r'))
        data=dict(np.load(fpath+'data.npz'))
        results=dict(np.load(fpath+'results.npz'))
        return cls(data=data,results=results,parameters=prm)

    def set_params(self,regen=0,test=0,**kwargs):
        """If some parameter 'name' in the template 'label' is changed by
        external arguments to the model of the form 'label_name', update it.
        If regen, then generate any matrix for which parameters havechanged
        In test mode, simply return which parameters would be changed by each kwarg."""
        if test or regen:
            changes={}
        for i,j in kwargs.iteritems():
            if i in self.dft_prm:
                if hasattr(self.dft_prm[i],'keys'):
                    self.prm.setdefault(i,{})
                    for k in j:
                        self.prm[i][k]=j[k]
                else:
                    self.prm[i]=j
                continue
            try:
                part=i.split('_')
                dico=self.prm
                p=part.pop(0)
                path=[]
                while part:
                    if p in dico:
                        path.append(p)
                        dico=dico[p]
                        p=part.pop(0)
                    else:
                        p+='_'+part.pop(0)
                if p in dico or kwargs.get('force',1) and p!=i:
                    if not p in dico:
                        print('set_params forcing:',i,p,j)
                    if test or (regen and j!=dico[p]):
                        changes[i]=path
                    if not test:
                        dico[p]=j
            except TypeError:
                pass
        if test:
            return changes
        if regen:
            #Only recreate whatever is changed by the imposed parameters
            print changes
            self.generate([i[0] for i in changes.values() ] )
        return self.prm


    def export_params(self,parameters=None):
        """Flatten parameter structure to the format 'label_name'. """
        if parameters is None:
            parameters=self.prm
        final={}
        def recurs(path,dico):
            if hasattr(dico,'keys'):
                for i,j in dico.iteritems():
                    recurs('{}_{}'.format(path,i),j )
            else:
                final[path]=dico
        for label in parameters:
            recurs(label,self.prm[label])
        return final



    def generate(self,labels=None):
        """Generate all the model coefficients based on parameters in self.prm. If labels, generate only those coefficients."""
        prm=self.prm
        N=prm['species']
        data=self.data
        dependencies={'growth':['environment','trophic','size'],'mortality':['environment'],'trophic_scale':['trophic','size'],
                      'trophic':['size'],
                      'competition_scale':['competition','size'], 'competition':['size']}

        def keep(name):
            # Must keep existing data with this name
            if labels is None:
                return 0
            if  name in labels or  set(labels).intersection(dependencies.get(name,[])) :
                return 0
            return 1

        def unwrap(dic,as_dict=1):
            """When nested dictionaries in prm (e.g. 'niche':{'width:{...}}, de-nest them into {'niche_width':...}"""
            if not hasattr(dic,'keys'):
                return []
            lst=[]
            for i in dic:
                k=dic[i]
                if hasattr(k,'keys'):
                    lst+=[(i,k)]+[('{}_{}'.format(i,j),l) for j,l in  unwrap(k,as_dict=0)]
            if as_dict:
                return dict(lst)
            return lst

        for name,dprm in unwrap(prm).iteritems():
            if keep(name):
                continue
            if not isinstance(dprm,dict) or not 'distribution' in dprm:
                continue
            if 'shape' in dprm:
                shape=dprm['shape']
            else:
                shape=N
            if not hasattr(shape,'__iter__'):
                shape=[shape]
            shape=[prm[sh] if isinstance(sh,basestring) else sh for sh in shape]

            dist=dprm['distribution']
            if dist=='uniform':
                res=np.random.uniform(dprm['range'][0],dprm['range'][1],shape )
            elif dist=='normal':
                if dprm['std']>0:
                    res=np.random.normal(dprm['mean'],dprm['std'],shape )
                else:
                    res=np.ones(shape)*dprm['mean']
            elif dist=='power':
                xmin,xmax=dprm['range']
                res= powlaw(dprm['exponent'],xmin,xmax,shape)
            elif dist=='noise':
                # Generate noisy landscape
                res=generate_noise(shape=shape,**{i:j for i,j in dprm.items() if i!='shape'} )
                rge=dprm['range']
                res=rge[0]+(res-np.min(res))*(rge[1]-rge[0])/(np.max(res)-np.min(res)) #Rescale generated noise

            if 'diagonal' in dprm:
                np.fill_diagonal(res,dprm['diagonal'])
            if 'sign' in dprm:
                sign=dprm['sign']
                if sign>0:
                    res=np.abs(res)
                elif sign<0:
                    res=-np.abs(res)
            if dprm.get('type','constant')=='variable':
                # Store variables in self.results
                self.results[name]=np.array([res])
            else:
                data[name]=res

            if name=='size':
                data['size']=np.sort(data['size']) #Order species by size (for convenience of plotting)


        trait=data['size']

        #=== Generate the energy structure (growth and trophic) ===

        # Generate trophic interactions with a niche model, using species body size as niche
        # Each predator is assigned an "eating range" (range of sizes it can eat) based on its own size

        growth = np.zeros(N)
        if not keep('trophic') or not keep('growth'):
            dprm = prm['trophic']
            dist = np.add.outer(-trait, trait).astype('float')  # Body size difference between predator and prey

            if dprm['ON']:
                if 'trophic' in data:
                    mat=data['trophic']
                else:
                    mat=np.ones((N,N) )

                if not keep('trophic'):
                    # Get center and width of eating range
                    center,width=dprm.get('distance',1),dprm.get('width',1)
                    center=center+width
                    range_exp=dprm.get('range_exp',0)
                    if not range_exp is 0:
                        oldcenter=center
                        center= center*(trait**range_exp).reshape( trait.shape+(1,) )
                        sig=dprm.get('randomness',0)
                        if sig>0:
                            center+= np.random.normal(0,sig, size=center.shape)
                        width=width*center/oldcenter
                        width[1:]=width[1:]+ np.min(np.abs(dist + center), axis=1).reshape(center.shape)[1:]
                              #np.min( [ ,  dw],axis=0)

                    # Set interactions to zero if the prey does not fall in the eating range
                    mat[dist > -center + width] = 0
                    mat[dist<-center-width ]=0
                    np.fill_diagonal(mat,0)
                    first=np.argmax(mat,axis=1)
                    notfirst=np.ones(mat.shape)
                    notfirst[(range(mat.shape[0]),first)]=0
                    nbpred=np.sum( np.sum(mat,axis=1)>0 )
                    mat[np.logical_and(np.random.random(mat.shape )>dprm.get('connectance',1),notfirst) ]=0
                    assert nbpred==np.sum( np.sum(mat,axis=1)>0 )
                    data['trophic']=mat
            else:
                mat=np.zeros((N,N))
            # Add autotrophic growth only to basal species
            growth[np.sum(mat,axis=1)==0 ]=1  #Species with no preys are autotrophs
            growth[trait==np.min(trait)]=1  #The smallest species is always an autotroph

        if not keep('trophic_scale'):
            data['trophic_scale']=prm['trophic']['scale']*trait**prm['trophic']['trait_exp'] #Spatial range scales with trait

        #=== Generate dispersal
        if not keep('dispersal'):
            data['dispersal']=prm['dispersal']['mean']*trait**prm['dispersal']['trait_exp']

        #=== Generate competition
        if not keep('competition_scale'):
            data['competition_scale']= prm['competition']['scale']*(trait**prm['competition']['trait_exp'])
        # print data['competition']

        #=== Generate growth and mortality
        env=data['environment']
        pos,wid,mortality=data['envniche_pos'],data['envniche_width'],data['mortality']
        pos,wid,growth,mortality= [z.reshape(z.shape+tuple(1 for i in range(len(env.shape))) ) for z in (pos,wid,growth,mortality) ]

        #Fitness with respect to abiotic environment, between 0 and 1, controls both growth and mortality
        abioticfit=np.exp(-(pos-env.reshape((1,)+env.shape))**2 /(2*wid**2)  )
        if not keep('growth'):
            data['growth']=growth=growth*abioticfit
        if not keep('mortality'):
            data['mortality']=mortality*(1-abioticfit)


        #=== Initial condition
        if not keep('n'):
            res=1
            if prm['n'].get('spatial',None) == 'spotty':
                dist = prm['n']['distribution_prm']
                width=dist.get('width',None)
                spots = range(dist['number'])
                axes = dist['axes']
                shape=self.results['n'][-1].shape
                othaxes = [a for a in range(len(shape)) if not a in axes]
                res = np.zeros(shape)
                for othidx in itertools.product(*[range(shape[oth]) for oth in othaxes]):
                    idx = [slice(None) for a in range(len(shape))]
                    for o, v in zip(othaxes, othidx):
                        idx[o] = v
                    idx = tuple(idx)
                    candidates = zip(*np.where(growth[idx]>.1))
                    candidates=[candidates[c] for  c in np.random.choice(range(len(candidates)),size=len(spots),replace=0)]
                    for c in candidates:
                        sidx=list(idx)
                        for o,v in zip(axes,c):
                            sidx[o]=v
                        sidx=tuple(sidx)
                        res[sidx] = 1
                    if not width is None:
                        kernel=np.ones([int(width*shape[a] ) for a in axes ])
                        if dist.get('shape',None)=='blot':
                            ksh=np.array(kernel.shape,dtype='int')
                            mn,sd=ksh/2.,ksh/4.
                            mn,sd=[x.reshape(x.shape+tuple(1 for i in range(len(ksh))) ) for x in (mn,sd) ]
                            kernel*=np.exp( - np.sum((np.indices(ksh) - mn)**2/(2*(sd)**2),axis=0) )
                            kernel[kernel<.2]=0
                            kernel/=np.max(kernel)
                        from scipy.signal import convolve
                        res[idx]=convolve(res[idx], kernel, mode='same')
            self.results['n'][-1]*=res


    def get_dx(self,t,x,calc_fluxes=False,eps=10**-6):
        """Get time derivatives  1/N dN/dt (used in differential equation solver below)

        Options:
            calc_fluxes: If True, do not only return the total derivative, but also the individual contribution
                of each term in the equation (dispersal, interactions...) to see which ones control the dynamics."""
        # x=x.copy()
        dx=np.zeros(x.shape)
        dxconst=np.zeros(x.shape)

        #Default convolution weights for dispersal when multiscaling switched off
        weights = np.zeros((3, 3), dtype='float')
        weights[0,1]=weights[2,1]=1.
        weights[1,0]=weights[1,2]=1.
        weights[1,1]=-np.sum(weights)

        data=self.data
        prm=self.prm
        mortality,growth,trophic,disp=data['mortality'],data['growth'],data['trophic'],data['dispersal']
        N=prm['species']
        death=prm.get('death',10**-15)
        dead=np.where(x<death)
        # x=np.clip(x,death,None)
        # x[dead]=0

        if calc_fluxes:
            typfluxes=['Trophic +', 'Trophic -', 'Source', 'Sink', 'Competition', 'Linear']
            fluxes=np.zeros((6,)+x.shape)

        rge=data['trophic_scale']
        for i in range(N):
            #Dispersal
            if prm['dispersal']['multiscale'] and 'dispersal_scale' in data:
                drge = data['dispersal_scale']
                xx = ndimage.gaussian_filter(x[i], sigma=drge[i], mode='wrap')
            else:
                xx=ndimage.convolve(x[i],weights,mode='wrap')
            dxdisp=disp[i]*xx
            dxconst[i]+=dxdisp
            if calc_fluxes:
                fluxes[2,i]+=np.clip(dxdisp,None,0)
                fluxes[3,i]+=np.clip(dxdisp,0,None)

            if not prm['trophic']['ON']:
                continue
            #Predation
            prey=np.where(trophic[i]!=0)[0]
            if not len(prey):
                continue
            if prm['trophic']['multiscale']:
                xx = ndimage.gaussian_filter(x[i], sigma=rge[i], mode='wrap')
                xps =[ ndimage.gaussian_filter(x[p], sigma=rge[i], mode='wrap')
                    for p in prey   ]
            else:
                xx= x[i]
                xps=[x[p] for p in prey   ]
            # pred = np.array([xx * xp *prm['trophic']['mean'] for xp  in xps])
            pred = np.array([xp *prm['trophic']['mean'] for xp  in xps])

            # dxprey=-pred*x[prey]/np.maximum(xps,eps)
            dxprey=-xx.reshape((1,)+xx.shape)
            dx[prey]+=dxprey
            # dxpred=np.sum(pred,axis=0)*x[i]/np.maximum(xx,eps) *prm['trophic']['efficiency'] #(1-x[i])
            dxpred=np.sum(pred,axis=0) *prm['trophic']['efficiency'] #(1-x[i])
            dx[i]+= dxpred
            if calc_fluxes:
                fluxes[0,i]+=np.abs(dxpred)
                fluxes[1,prey]+=-np.abs(dxprey)

        #Competition
        if prm['competition']['ON']:
            comp=data['competition']
            xc = x
            if prm['competition']['multiscale']:
                rge = data['competition_scale']
                xc = np.array([ndimage.gaussian_filter(x[i], sigma=rge[i], mode='wrap') for i in range(N)])
            dxcomp= np.tensordot(comp,xc, axes=(1,0))
        else:
            dxcomp=0

        dxlin=growth-mortality
        dx += dxlin-dxcomp

        dx[dead]=np.clip(dx[dead],0,None)
        if calc_fluxes:
            fluxes[4]+=-np.abs(dxcomp)
            fluxes[5]+=dxlin
            return dx,typfluxes,fluxes
        if t==0 and 0:
            plt.figure()
            plt.imshow(x[0]*dxlin[0])
        dx=dxconst + (x-death)*dx
        return dx

    def get_kernel_FT(self,a,k,kernel='gauss'):
        """Fourier Transform of the spatial kernel"""
        lx, ly = k.shape
        if a<=0.1:
            return np.ones(k.shape)
        if kernel=='gauss':
            val= np.exp(-(k**2)*(a**2/2))#*np.pi*(a**2*2) /(  np.pi * drge[i] ** 2 * 2
            return val
        elif kernel=='exp':
            return 2*a/(a**2+np.abs(k)**2)

    def get_dx_FT(self, t, x, calc_fluxes=0,resconv=128):
        """Get time derivatives in Fourier space (used in differential equation solver below) """
        dx = np.zeros(x.shape,dtype='complex')
        dxdisp = np.zeros(x.shape,dtype='complex')
        size=np.prod(dx[0].shape).astype('float')

        data = self.data_FT
        prm = self.prm
        mortality, growth, trophic, disp = data['mortality'], data['growth'], data['trophic'], data['dispersal']
        k,kf=data['knormFourier'],data['kFourier']
        N = prm['species']

        if calc_fluxes:
            typfluxes = ['Trophic +', 'Trophic -', 'Dispersal', 'Competition', 'Linear']
            fluxes = np.zeros((5,) + x.shape)

        for i in range(N):
            # Dispersal
            dxdisp[i] = - k**2 * disp[i] * x[i]
            if calc_fluxes:
                fluxes[2, i] += np.abs(dxdisp)

            if not prm['trophic']['ON']:
                continue
            #Predation
            prey = np.where(trophic[i] != 0)[0]
            if not len(prey):
                continue

            Aij=np.tensordot(trophic[i][prey] , self.get_kernel_FT(self.data['trophic_scale'][i],k,
                                                                prm['trophic']['kernel'] ),axes=0)

            dxprey = - x[i] * Aij
            dx[prey] += dxprey
            dxpred = np.sum(Aij*x[prey], axis=0 ) * prm['trophic']['efficiency']
            dx[i] += dxpred
            if calc_fluxes:
                fluxes[0, i] += np.abs(dxpred)
                fluxes[1, prey] += np.abs(dxprey)

        # Competition
        if prm['competition']['ON']:
            comp = data['competition']
            dxcomp = np.tensordot(comp, x*[ self.get_kernel_FT(self.data['competition_scale'][i],k,
                                                    prm['competition']['kernel']) for i in range(N)], axes=(1, 0))
        else:
            dxcomp=0
        dxlin = growth - mortality
        dx += dxlin - dxcomp

        if calc_fluxes:
            fluxes[3] += np.abs(dxcomp)
            fluxes[4] += np.abs(dxlin)
            return dx, typfluxes, fluxes

        def window(z,resconv):
            # Only convolve with a limited window; typically remove at least the edges
            if resconv is None:
                return z
            lx,ly=z.shape
            cx,cy=lx/2+1,ly/2+1
            z2= z[max(cx-resconv,0):cx+resconv, max(cy-resconv,0):cy+resconv]
            return z2
        return np.array([ssignal.fftconvolve( dx[i]/ size,window(x[i],min(resconv,min(dx.shape[1:])//2-1 ))
                                                                 ,'same') for i in range(N)]) + dxdisp


    def prep_FT(self,x,**kwargs):
        """Prepare data for Fourier-transformed dynamical equations (if used)"""

        use_Fourier = kwargs.get('use_Fourier',1)
        if not use_Fourier:
            return x
        data_FT={}
        data=self.data
        landscape=data['environment']
        lx,ly=landscape.shape
        kx,ky=fft.fftfreq(lx)*2*np.pi, fft.fftfreq(ly)*2*np.pi
        for i in ('mortality','growth','trophic','competition','dispersal'):
            if i in data and data[i].shape[1:]==landscape.shape:
                data_FT[i]= reorder_FT( fft.fft2(data[i]), (kx,ky))
                assert checkpos_FT(data_FT[i])<10**-10
            else:
                data_FT[i]=data[i]
        if use_Fourier=='reduced':
            kx=kx[:lx/2]
            ky=ky[:ly/2]
        data_FT['kxFourier'],data_FT['kyFourier']= kx*lx,ky*ly
        data_FT['knormFourier']= reorder_FT(np.sqrt(np.add.outer(kx**2,ky**2 )),(kx,ky))
        data_FT['kFourier']= reorder_FT(np.sqrt(np.add.outer((kx*lx)**2,(ky*ly)**2 )),(kx,ky))
        self.data_FT=data_FT
        x2=fft.fft2(x)
        x2= reorder_FT(x2,(kx,ky))
        assert checkpos_FT(x2) <10**-10
        return x2

    def save_results(self,t,x,use_Fourier=0,print_msg=0,death=0):
        if use_Fourier:
            kx, ky = self.data_FT['kxFourier'], self.data_FT['kyFourier']
            if use_Fourier == 'reduced':
                lx, ly = self.data['environment'].shape
                x2 = np.zeros(x.shape)
                x2[:lx / 2, :ly / 2] = x
                x2[lx / 2:, :ly / 2] - np.conj(x)
                x2[:lx / 2, ly / 2:] - np.conj(x)
                x2[lx / 2:, ly / 2:] = x
            else:
                x2 = x
            realx = fft.ifft2(reorder_FT(x2, (kx, ky))).real
        else:
            realx = x
        realx[realx <= death] = 0
        if print_msg:
            print('   Btot {:.2g}'.format(np.sum(realx)))
        self.results['n'] = np.concatenate([self.results['n'], [realx]])
        self.results['t'] = np.concatenate([self.results['t'], [t]])

    def integrate(self,integrator,ti,tf,x,use_Fourier=0,**kwargs):
        tries=0
        if np.sum(x)<kwargs.get('death',10**-10):
            #If everyone is dead, just skip
            return  x, True, '',None
        deltat=kwargs.get('deltat',None)
        if deltat is None:
            deltat= tf - ti
        t=lastfail=ti
        result,success,error=None,False,None
        while (not success) or t<tf:
            if tries>kwargs.get("MAX_TRIES",500):
                error="ERROR: MAX TRIES reached in LandscapeModel.integrate"
                result=x
                t=tf
                break
            if tries>0 and deltat>10**-15:
                lastfail=t
                if deltat/2 >= kwargs.get('mindeltat',0):
                    deltat/=2.
                    if kwargs.get('print_msg',1):
                        print 'Failure, trying tsample', deltat
                else:
                    error = 'ERROR: step too small for integrator to converge'
                    break
            integrator.set_initial_value(x.ravel(), t)  # .set_f_params(reff)
            result=integrator.integrate(t+max(10**-5,min(deltat,tf-t)) )
            success=integrator.successful()
            if success:
                t = integrator.t
                x = result
                if t<tf and not tries and t>deltat*100 and t-lastfail>deltat*10 :
                    deltat*=2
            tries+=1
        if np.isnan(result).any():
            error = 'ERROR: NaN'
        return result,success, error,ifelse(deltat<(tf-ti),deltat,None)

    def evol(self,tmax=10,nsample=10,dt=.1,keep='all',print_msg=1,method='scipy',use_Fourier=0, **kwargs):
        """Time evolution of the model"""

        if not method in self.methods:
            raise Exception('ERROR: Method "{}" not found'.format(method))

        if kwargs.get('reseed',1) or not self.data:
            # Create new matrices and initial conditions
            self.generate()
            self.results['t']=np.zeros(1)
        else:
            init = kwargs.get('init',0)
            if init=='extend':
                # Restart from end of last simulation
                pass
            elif init=='restart':
                # Restart from start of last simulation
                for label in self.results:
                    self.results[label]=self.results[label][:1]
            else:
                # Only create new initial conditions
                print "Regenerating initial conditions"
                self.generate(self.results.keys() )
                self.results['t']=np.zeros(1)

        x0=x=self.results['n'][-1].copy()

        death = self.prm.get('death', 10 ** -15)

        if use_Fourier:
            x=self.prep_FT(x,**kwargs)
            lx,ly=x.shape[-2:]
            cx,cy=lx/2,ly/2

            def get_dx(t,x):
                x=x.reshape(x0.shape)
                return self.get_dx_FT(t,x ).ravel()
            integ='zvode'
        else:
            if 'noise' in self.data:
                integ = 'lsoda'
            else:
                integ = 'dop853'
            def get_dx(t,x):
                x=np.clip(x,0,None).reshape(x0.shape)
                return (self.get_dx(t,x)).ravel()

        t,deltat=0,None
        if kwargs.get('samplescale','log'):
            tsamples=list(np.logspace(np.log10(dt),np.log10(tmax),nsample ))
        else:
            tsamples = list(np.linspace(dt,tmax, nsample))

        if method=='scipy':
            integrator = scint.ode(get_dx).set_integrator(integ, nsteps=500000)
            for ts in tsamples:
                x, success, error,deltat = self.integrate(integrator,t,ts, x, use_Fourier=use_Fourier,print_msg=print_msg,
                                                          deltat=deltat,**kwargs)
                if error:
                    print error
                    return 0
                if not success:
                    print 'WARNING: scipy integrator failed, switching to Euler.'
                    method='Euler'
                    break
                t=ts
                if print_msg:
                    print('Time {}'.format(ts) )
                if keep=='all' or ts+dt>=tmax:
                    xx=x.reshape(x0.shape)
                    if np.max(np.abs(xx-self.results['n'][-1]))<10**-5:
                        print 'WARNING: EQUILIBRIUM REACHED'
                        break
                    self.save_results(t,xx,use_Fourier=use_Fourier,print_msg=print_msg,death=death)
        if method=='Euler':
            while t<tmax:
                if  use_Fourier:
                    dx=self.get_dx_FT(t,x)
                    x=setpos_FT(x+dt*dx)
                    # x=x+dt*dx
                    x[:,cx,cy]=np.clip(x[:,cx,cy],0,None)
                else:
                    dx=self.get_dx(t,x)
                    x+=dt* dx
                    x[x<10**-15]=0
                t+=dt
                if  t+dt > tsamples[0]:
                    tsamp=tsamples.pop(0)
                    if not tsamples:
                        tsamples.append(tmax)
                    if print_msg:
                        print('Time {}'.format(t) )
                    if keep=='all' or t+dt>tmax:
                        self.save_results(tsamp, x, use_Fourier=use_Fourier, print_msg=print_msg,death=death)

        return 1



#===========================================================================================================

# LOOPER

#===========================================================================================================

def mp_run(array):
    """Parallel processing utility function."""
    path,Model,prm,results,tmax,tsample,keep,replica=array
    model=Model(results=results,**prm)
    model.PARALLEL_LOCK=True
    model.evol(tmax=tmax,tsample=tsample,keep=keep)
    model.save(path,suffix=replica )
    model.module['prey'].tree=None

class Looper(object):
    """Allows to recursively loop over parameter range, running a Model then saving results in a tree of folders."""
    MAX_TRIALS=1
    def __init__(self,Model):
        self.Model=Model
        self.model={}

    def loop_core(self,tsample=50.,tmax=None, path='.',**kwargs):
        """Main function (internal use)."""
        import pandas as pd
        path=Path(path)
        path.norm()
        path.mkdir(rmdir=kwargs.get('clean_directory',False))
        table=[]

        dic={
            'tsample':tsample,
            'tmax':tmax,
            'path':path,
            'sys':kwargs.get('sys',None),
            }

        reseed=kwargs.get('reseed',1)
        if self.model and not reseed:
            model=deepcopy(self.model[dic['sys']])
            model.set_params(regen=1,**kwargs)
        else:
            model=self.Model(**kwargs)
        for i in kwargs:
            if i in model.export_params():
                dic[i]=kwargs[i]

        if not reseed:
            self.model[dic['sys']]=model #transmit the model further down the loop instead of recreating it each time

        if 'replicas' in kwargs:
            '''Parallel processing'''
            from multiprocessing import Pool,freeze_support
            freeze_support()
            systems=range(kwargs['replicas'])

            array=[(path.copy(),self.Model,model.export_params(),model.results,tmax,tsample,kwargs.get('keep','endpoint'),sys)
                     for sys in systems  ]

            pool=Pool()
            pool.map(mp_run,array )
            pool.close()
            pool.join()
            pool.terminate()
        else:
            kwargs.setdefault('converge',1)
            kwargs.setdefault('keep','endpoint')
            success=model.evol(tmax=tmax,tsample=tsample,path=path,**kwargs)
            if not success and not 'model' in kwargs:
                kwargs['loop_trials']=kwargs.get('loop_trials',0)+1

                if kwargs['loop_trials']<self.MAX_TRIALS:
                    print('### model diverged, trying again! ###')
                    return self.loop_core(tsample=tsample,tmax=tmax, path=path,**kwargs)
                else:
                    print('### model diverged {} times, setting results to 0! ###'.format(kwargs['loop_trials']))
                    for var in model.results:
                        model.results[var][:]=0
            if 'loop_trials' in kwargs:
                del kwargs['loop_trials']
            model.save(path )

        table.append(dic)
        pd.DataFrame(table).to_csv(path+'files.csv')
        return table

    def loop(self,axes=None,path='.',*args,**kwargs):
        """External loop function: give it a list of parameters and parameter values to iterate over as "axes",
        and a path to save the results in. """
        respath=Path(path)

        if axes is None or len(axes)==0:
            return []
        axis=axes[0]
        axes=axes[1:]


        key=axis[0]
        if 'check_key' in kwargs and not kwargs['check_key']==False:
            if kwargs['check_key'] in (True,1):
                kwargs['check_key'] = self.Model(**kwargs).export_params()
            to_check=to_iterable(key)
            for c in to_check:
                if not c in kwargs['check_key']:
                    raise Exception('Key "{}" not found among possible parameters.'.format(c))

        table =[]
        labels=[]
        refkwargs=kwargs
        for x in axis[1]:
            kwargs={}
            kwargs.update(refkwargs) #Prevents changes to kwargs from carrying over

            if hasattr(key,'__iter__'):
                #Case where multiple factors must change jointly (CANNOT BE SYS OR INIT)
                for i,j in zip(key,x):
                    kwargs[i]=j
                folder=Path('{}-{}_etc'.format(key[0],x[0]) )
                base=str(folder)
                inti=2
                while str(folder) in labels:
                    folder=Path(base+'{}'.format(inti))
                    inti+=1
                labels.append(str(folder))

            else:
                kwargs[key]=x
                folder=Path('{}-{}'.format(key,x) )

            if axes or not key in ('init','replica','sys'):
                #Do not print the innermost level of recursion if it is either
                #initial conditions, or replica (for stochastic systems)
                if kwargs.get('print_msg',1):
                    print(folder)

            kwargs.setdefault('init','restart')
            if axes:
                table+=self.loop(*args,axes=axes,path=respath+folder, **kwargs)
            else:
                fpath=Path(respath + folder)
                fpath.norm(),fpath.mkdir()
                if not kwargs.get('rerun',0) and 'files.csv' in os.listdir(fpath):
                    table+=[x.to_dict() for idx,x in pd.read_csv(fpath+'files.csv',index_col=None).iterrows()]
                else:
                    table+=self.loop_core(*args,path=fpath, **kwargs)
        return table

