from utilities import *

class LandscapeModel():
    """Main Model object: given parameters, will generate matrices, run dynamics, load and save."""

    #Default parameters
    dft_prm={
        'species':30,

        #Landscape size
        'landx':100,
        'landy':100,

        #Abundance of each species at each pixel
        'n': {
            'type':'variable',
            'shape': ('species','landx','landy'),
            'range': (0., 1.),
            'distribution': 'uniform',
        },

        #Value of the environment at each pixel
        'environment':{
            'shape':('landx','landy'),
            'range':(0.,100.),
            'distribution':'noise',
            'spectralexp':-1,
        },

        #Species log-body size (trait used below to parametrize all interactions, trophic and nontrophic)
        'size': {
            'range': (0, 100),
            'distribution': 'uniform',
        },

        #Environmental niche
            #Center
        'envniche_pos':{
            'range':(20,80),
            'distribution':'uniform',
        },
            #Width
        'envniche_width':{
            'mean':20., 'std':5,
            'distribution':'normal',
        },

        #Trophic interactions
        'trophic':{
            'mean':0.5,
            'efficiency':0.9,
            'distance':.5,
            'width':.5,
            'rangeexp':1, #Exponent for how eating range (size range of edible prey) changes with trait value
            'multiscale':0, #Switch ON/OFF multiscale
            'traitexp': 1, #Exponent for how spatial range of interaction changes with trait value
        },

        #Dispersal
        'dispersal':{
            'mean':0.01,
            'multiscale':0,  #Switch ON/OFF multiscale
            'traitexp': 1, #Exponent for how spatial range of interaction changes with trait value

        },

        #Competition
        'competition': {
            'mean': 0.1,
            'multiscale': 0, #Switch ON/OFF multiscale
            'traitexp':1, #Exponent for how spatial range of interaction changes with trait value
        },

        #Mortality
        'mortality': {
            'distribution': 'uniform',
            'range': (0.001, 0.01),
        },

    }


    def __init__(self,**kwargs):
        """When the model object is created, store passed information."""
        self.data=deepcopy(kwargs.pop('data',{}))
        self.results=deepcopy(kwargs.pop('results',{}))
        self.prm=deepcopy(kwargs.pop('parameters',LandscapeModel.dft_prm))
        self.set_params(**kwargs)


    def save(self,path):
        """Save model parameters in prm.dat, dynamical variables in results.npz and other matrices in data.npz"""
        fpath=Path(path)
        fpath.mkdir()
        dumps(open(fpath+'prm.dat','w'),self.prm)
        np.savez(fpath+'results',**self.results)
        np.savez(fpath+'data',**self.data)


    @classmethod
    def load(klass,path):
        """Load variables and other parameters from save files in path."""
        fpath=Path(path)
        prm=loads(open(fpath+'prm.dat','r'))
        results=dict(np.load(fpath+'results.npz'))
        data=dict(np.load(fpath+'data.npz'))
        return klass(data=data,results=results,parameters=prm)

    def set_params(self,**kwargs):
        """If some parameter 'name' in the template 'label' is changed by
        external arguments to the model of the form 'label_name', update it."""
        for i,j in kwargs.iteritems():
            if i in self.dft_prm:
                self.prm.setdefault(i,{})
                for k in j:
                    self.prm[i][k]=j[k]
                continue
            try:
                part=i.split('_')
                dico=self.prm
                p=part.pop(0)
                while part:
                    if p in dico:
                        dico=dico[p]
                        p=part.pop(0)
                    else:
                        p+='_'+part.pop(0)
                if p in dico or kwargs.get('force',1) and p!=i:
                    if not p in dico:
                        print('set_params forcing:',i,p,j)
                    dico[p]=j
            except:
                pass
        return self.prm

    def generate(self):
        """Generate all the matrices based on parameters self.prm."""
        prm=self.prm
        data=self.data
        N=prm['species']

        #=== Generate random distributions ===

        for name in prm:
            dprm=prm[name]
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
                res=np.random.normal(dprm['mean'],dprm['std'],shape )
            elif dist=='noise':
                # Generate noisy landscape
                res=generate_noise(shape=shape,**{i:j for i,j in dprm.items() if i!='shape'} )
                rge=dprm['range']
                res=rge[0]+res*(rge[1]-rge[0])/(np.max(res)/np.min(res)) #Rescale generated noise

            if dprm.get('type','constant')=='variable':
                # Store variables in self.results
                self.results[name]=np.array([res])
            else:
                data[name]=res
        data['size']=np.sort(data['size']) #Order species by size (for convenience of plotting)


        #=== Generate the energy structure (growth and trophic) ===

        # Generate trophic interactions with a niche model, using species body size as niche
        # Each predator is assigned an "eating range" (range of sizes it can eat) based on its own size
        trait=data['size']
        dprm=prm['trophic']
        dist=np.add.outer(-trait,trait).astype('float') # Body size difference between predator and prey
        mat=np.ones((N,N) )

        # Get center and width of eating range
        center,width=dprm.get('distance',1),dprm.get('range')
        if dprm.get('rangeexp',0)!=0:
            oldcenter=center
            center= center*trait**rangeexp.reshape( trait.shape+(1,) )
            width=width*center/oldcenter+np.min(np.abs(dist),axis=1)

        # Set interactions to zero if the prey does not fall in the eating range
        mat[dist > -center + width] = 0
        mat[dist<-center-width ]=0
        np.fill_diagonal(mat,0)
        data['trophic']=mat

        # Add heterotrophic growth only to basal species
        growth=np.zeros(N)
        growth[np.sum(mat,axis=1)<1 ]=1  #Species with no preys are heterotrophs
        growth[trait==np.min(trait)]=1  #The smallest species is also a heterotroph
        data['trophic_range']=trait**prm['trophic']['traitexp'] #Spatial range scales with trait

        #=== Generate dispersal
        data['dispersal']=trait**prm['dispersal']['traitexp']

        #=== Generate competition
        data['competition_range']=trait**prm['dispersal']['traitexp']
        data['competition']=np.eye(N)

        #=== Generate growth and mortality
        env=data['environment']
        pos,wid,mortality=data['envniche_pos'],data['envniche_width'],data['mortality']
        pos,wid,growth,mortality= [z.reshape(z.shape+tuple(1 for i in range(len(env.shape))) ) for z in (pos,wid,growth,mortality) ]

        #Fitness with respect to abiotic environment, between 0 and 1, controls both growth and mortality
        abioticfit=np.exp(-(pos-env.reshape((1,)+env.shape))**2 /(2*wid**2)  )
        data['growth']=growth*abioticfit
        data['mortality']=mortality*(1-abioticfit)

    def get_dx(self,t,x,calc_fluxes=False):
        """Get time derivatives (used in differential equation solver below)

        Options:
            calc_fluxes: If True, do not only return the total derivative, but also the individual contribution
                of each term in the equation (dispersal, interactions...) to see which ones control the dynamics."""

        dx=np.zeros(x.shape)


        #Default convolution weights when multiscaling switched off
        weights = np.ones((3, 3), dtype='float')
        weights /= np.sum(weights)

        data=self.data
        prm=self.prm
        mortality,growth,inter,disp=data['mortality'],data['growth'],data['trophic'],data['dispersal']
        N=prm['species']
        death=prm.get('death',10**-15)
        dead=np.where(x<death)
        x=np.clip(x,death,None)

        if calc_fluxes:
            typfluxes=['Trophic +', 'Trophic -', 'Dispersal', 'Competition', 'Linear']
            fluxes=np.zeros((5,)+x.shape)

        #Predation
        rge=data['trophic_range']
        for i in range(N):
            prey=np.where(inter[i]!=0)[0]
            if not len(prey):
                continue
            if prm['trophic']['multiscale']:
                xx= ndimage.gaussian_filter(x[i], sigma=rge[i] )
                xps=[ndimage.gaussian_filter(x[p],sigma=rge[i] ) for p in prey   ]
            else:
                xx= x[i]
                xps=[x[p] for p in prey   ]
            pred = np.array([xx * xp *prm['trophic']['mean'] for xp  in xps])

            dxprey=-pred*x[prey]/xps
            dx[prey]+=dxprey
            dxpred=np.sum(pred,axis=0)*x[i]/xx*(1-x[i]) *prm['trophic']['efficiency']
            dx[i]+= dxpred
            if calc_fluxes:
                fluxes[0,i]+=np.abs(dxpred)
                fluxes[1,prey]+=np.abs(dxprey)

            #Dispersal
            if prm['dispersal']['multiscale']:
                xx= ndimage.gaussian_filter(x[i], sigma=disp[i] )
            else:
                xx=ndimage.convolve(x[i],weights,mode='wrap')
            # else:
            dxdisp=prm['dispersal']['mean']*(xx-x[i])
            # code_debugger()
            dx[i]+=dxdisp

            if calc_fluxes:
                fluxes[2,i]+=np.abs(dxdisp)

        #Competition
        comp=data['competition']
        xc = x
        if prm['competition']['multiscale']:
            rge = data['competition_range']
            xc=ndimage.gaussian_filter(x, sigma=rge)
        dxcomp= x*np.tensordot(comp,xc, axes=(1,0))
        dxlin=growth-mortality
        dx += dxlin-dxcomp

        dx[dead]=np.clip(dx[dead],0,None)
        if calc_fluxes:
            fluxes[3]+=np.abs(dxcomp)
            fluxes[4]+=np.abs(dxlin)
            return dx,typfluxes,fluxes
        return dx

    def evol(self,tmax=5,tsample=.1,dt=.1,keep='all',print_msg=1,**kwargs):
        """Time evolution of the model"""

        #Create new matrices and initial conditions
        if kwargs.get('reseed',1):
            self.generate()

        x=self.results['n'][-1]

        t=0
        while t<tmax:

            dx=self.get_dx(t,x)
            t+=dt
            x*=(1+np.clip(dt* dx,-.999,None))
            x[x<10**-10]=0

            if t%tsample <dt or t+dt > tmax:
                if print_msg:
                    print('Time {}'.format(t) )
                if keep=='all' or t+dt>tmax:
                    self.results['n']=np.concatenate([self.results['n'],[x]])

        return 1



class Looper(object):
    """Allows to recursively loop over parameter range, running a Model then saving results in a tree of folders."""

    def __init__(self,Model):
        self.Model=Model

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


        if 'model' in kwargs:
            model=kwargs['model']
            model.set_params(**kwargs)
        else:
            model=self.Model(**kwargs)
        if not kwargs.get('reseed',1):
            kwargs['model']=model #transmit the model further down the loop

        # dic.update(model.export_params())

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
                        model.results[var].matrix[:]=0
            if 'loop_trials' in kwargs:
                del kwargs['loop_trials']
            model.save(path,overwrite=1 )

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


            if axes:
                table+=self.loop(*args,axes=axes,path=respath+folder, **kwargs)
            else:
                table+=self.loop_core(*args,path=respath+folder, **kwargs)
        return table

