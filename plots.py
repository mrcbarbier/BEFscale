from landscapemodel import *
import matplotlib.figure as mpfig


# ============= MEASURES =================

def basic_measures(model,dic):
    """Measure various quantities for a given model, store summary in dictionary dic."""
    Nf = model.results['n'][-1]
    growth=model.data['growth']

    B=np.sum(Nf,axis=0)  #Total Biomass per patch
    basal=np.where(np.max(growth,axis=(1,2))>0 )[0] #Basal species
    Bbasal=np.sum(Nf[basal ],axis=0)  #Total producer biomass per patch
    prod=np.sum(growth*Nf,axis=0)  #Total production per patch

    alive=(Nf > model.prm['death'])
    n=Nf/(B+10**-15) # Relative biomass
    D = np.sum(alive, axis=0)  # Diversity
    DS = -np.sum(n * np.log(np.clip(n, 10 ** -15, None)), axis=0)  # Shannon diversity
    dic.update({ 'biomass_tot':np.sum(B),'biomass_basal':np.sum(Bbasal),
                 'production_tot':np.sum(prod),
                 'alpha_div':np.mean(D), 'alpha_div_min':np.min(D), 'alpha_div_max':np.max(D),
                 'alpha_shannon':np.mean(DS), 'alpha_shannon_min':np.min(DS), 'alpha_shannon_max':np.max(DS),
                 'gamma_div':np.sum(np.sum(alive,axis=(1,2))>0 ) })


    # from scipy.signal import correlate2d
    # env=model.data['environment']
    # dic.update({'landscape_scale':correlate2d(env,env)})

    return B, Bbasal, D, DS

# ============= PLOTS =================

def detailed_plots(path,save=0,movie=0,**kwargs):
    """Plots with detailed information on each simulation outcome"""

    if save:
        if isinstance(save, basestring) or isinstance(save, Path):
            outpath = save
        else:
            outpath = path + Path('plots')

    files=pd.read_csv(path+'files.csv')
    table=[]
    for idx,f in files.iterrows():
        model=LandscapeModel.load(f['path'])
        print 'Plotting',f['path']
        dic = dict(model.export_params())
        dic['path'] = f['path']
        B, Bbasal, D, DS=basic_measures(model, dic)
        table.append(dic)
        figs = {}

        prm,data=model.prm,model.data
        death=prm['death']

        if save:
            fpath=Path(outpath+'_'.join([n for n in Path(f['path']).split() if not n in Path(outpath).split() ]))
            fpath.mkdir()

        nf=model.results['n'][-1]
        N=prm['species']
        size=model.data['size']

        #Species trajectories
        figs['N_traj']=plt.figure()
        nsum=np.sum(model.results['n']/death, axis=(2, 3))
        for i in zip(*nsum):
            plt.plot(model.results['t'],i)
        plt.yscale('log')
        plt.xscale('log')
        plt.xlabel('Time')
        plt.title('Total abundance per species')

        #Plots of species abundance per patch
        if save and movie:
            ts=range(model.results['n'].shape[0])
        else:
            ts=[0,-1]
        for t in ts:
            fig = plt.figure(figsize=np.array(mpfig.rcParams['figure.figsize']) * (3, 3))
            if t==0:
                figs['N_initial']=fig
                plt.suptitle('Initial abundance')
            panel=0
            for i in range(N):
                panel,ax=auto_subplot(panel,N)
                plt.imshow(model.results['n'][t][i]/death )#,vmin=0,vmax=1.2)
                plt.colorbar()
                plt.title('Species {}'.format(i))
            if save and movie:
                print 'Movie frame {}'.format(t)
                mpath = fpath + Path('movie')
                fig.savefig(mpath.mkdir() + 'fig{}.png'.format(t))
                if t != ts[0] and t != ts[-1]:
                    plt.close(fig)

        figs['N_final']=fig
        plt.suptitle('Final abundance')


        #Species environmental niches
        # figenv=plt.figure()
        # panel=0
        figs['Environment'], (a0, a1,a2) = plt.subplots(1, 3, gridspec_kw={'width_ratios': [6, 1,4]},
                                           figsize=np.array(mpfig.rcParams['figure.figsize']) * (2, 1))

        # panel, ax = auto_subplot(panel, 3,rows=1)
        #plt.colorbar(
        a0.imshow(data['environment'])#,ax=a0)
        a0.set_xticks([]),a0.set_yticks([])
        a0.set_title('Landscape')
        # plt.colorbar()
        # panel, ax = auto_subplot(panel, 3,rows=1)
        a1.hist(data['environment'].ravel(),bins=20,
                       orientation=u'horizontal',  )
        a1.set_xticks([])
        a1.set_title('Env.')
        #ax.set_aspect(4.)
        # panel, ax = auto_subplot(panel, 3,rows=1)

        c,w=data['envniche_pos'], data['envniche_width']
        a2.vlines(np.arange(N),c-w,c+w )
        a2.set_title('Species abiotic niche')
        plt.ylim(ymin=np.min(data['environment']),ymax=np.max(data['environment']))
        figs['Environment'].tight_layout()


        #Food web
        if model.prm['trophic']['ON']:
            Pij=model.data['trophic']
            trophic_height=naive_trophic_score(Pij,(np.sum(model.data['growth'],axis=(1,2))>0).astype('float') )+1
            figs['foodweb']=plt.figure(figsize=np.array(mpfig.rcParams['figure.figsize'])*(2,2) )
            draw_network(Pij,#xpos=trophic_height/(size/np.min(size)+1)*np.random.normal(1,.2,size=size.shape),
                         ypos=size,newfig=0,hold=1)
            plt.ylabel('Body size')
            plt.title('Food web')

        #Role of different interactions
        figs['dynamics']=plt.figure(figsize=np.array(mpfig.rcParams['figure.figsize'])*(3,3) )
        panel=0
        dx,typfluxes,fluxes=model.get_dlogx(None, nf, calc_fluxes=1)
        for i in range(N):
            panel,ax=auto_subplot(panel,N)
            for t, typ in enumerate(typfluxes):
                ff=fluxes[t,i]/np.clip(np.abs(dx[i]),10**-10,None)
                # code_debugger()
                ff=np.mean(ff[nf[i]>death] )  #Average over all patches
                # f=f.ravel()[np.argmax(nf[i])]  #Values at most abundant patch
                plt.bar(t,ff)
            plt.title('Species {}'.format(i))
        plt.legend(typfluxes)
        plt.suptitle('Dynamical effect of interactions')


        #Biomass pyramid
        figs['N_pyramid']=plt.figure()
        ys=np.sum(nf/death,axis=(1,2))
        survivors=(ys>1)
        plt.scatter(data['size'][survivors],ys[survivors])
        plt.xlabel('Size')
        plt.yscale('log')
        plt.xscale('log')
        plt.ylabel('Abundance')
        plt.title('Size spectrum')

        #SAR
        figs['SAR']=plt.figure()
        nstart=200
        startpos=(np.random.randint(prm['landx'],size=nstart),np.random.randint(prm['landy'],size=nstart) )
        dist=(np.add.outer(startpos[0],-startpos[0])**2+np.add.outer(startpos[1],-startpos[1])**2)
        SARS=[[],[]]
        for i in range(nstart):
            order=np.argsort(dist[i])
            xs=dist[i,order]
            ys=[nf[:,startpos[0][o],startpos[1][o]]>death for o in order]
            ys=np.cumsum(ys,axis=0)
            ys=np.sum(ys>0,axis=1)
            plt.plot(xs,ys,alpha=0.2)
            SARS[0]+=list(xs)
            SARS[1]+=list(ys)
            # code_debugger()
        SARS=np.array(SARS)
        xs,ys=SARS
        bins=np.logspace(np.min(np.log10(xs[xs>0])),np.max(np.log10(xs[xs>0])),20)
        ibins=np.digitize(xs,bins)
        nonempty=[i for i in range(len(bins)) if (ibins==i).any() ]
        sxs,sys=bins[nonempty],[np.mean(ys[ibins==i] )  for i in nonempty ]
        plt.plot(sxs,sys ,lw=3,c='k')
        dic['SAR']=np.array([sxs,sys])
        plt.xscale('log')
        plt.title('SAR')
        plt.xlabel('Area')
        plt.ylabel('Diversity')
        # dic.update({'min':strength})


        #BEF
        figs['BEF']=plt.figure()
        patches=pd.DataFrame([dict( zip(('B','Bbasal','D','DShannon'), x)  ) for x in zip(B.ravel(),Bbasal.ravel(),
                                                                                          D.ravel(),DS.ravel()) ]  )

        panel=0
        for name,val in [('Biomass','B'),('Basal biomass','Bbasal')]:
            panel,ax=auto_subplot(panel,4)
            # Group by species diversity
            patchesD=patches.groupby('D').median()
            plt.scatter(patches['D'],patches[val])
            plt.yscale('log')#,plt.xscale('log')
            plt.plot(patchesD.index,patchesD[val],lw=2,c='k' )
            plt.xlabel('Diversity')
            plt.ylabel(name)
            panel,ax=auto_subplot(panel,4)

            # Group by Shannon
            bins = np.linspace(patches['DShannon'].min(), patches['DShannon'].max(), 5)
            patchesDS=patches.groupby(np.digitize(patches['DShannon'], bins)).median()
            plt.scatter(patches['DShannon'],patches[val])
            plt.yscale('log')#, plt.xscale('log')
            plt.plot(patchesDS['DShannon'],patchesDS[val],lw=2,c='k' )
            plt.xlabel('Shannon diversity')
            plt.ylabel(name)

        #Checkerboard pattern
        figs['checkerboard']=plt.figure(figsize=np.array(mpfig.rcParams['figure.figsize']) * (3, 1))
        corr=np.corrcoef([nf[i,:,:].ravel() for i in range(nf.shape[0])])
        plt.subplot(131)
        plt.colorbar(plt.imshow(corr, cmap='seismic_r', vmin=-1, vmax=1),
                     ax=plt.gca())
        plt.title('Spatial correlation')
        plt.subplot(132)
        mat=-model.data['competition']
        if prm['trophic']['ON']:
            tmat=model.data['trophic']
            mat-=tmat.T-prm['trophic']['efficiency']*tmat
        plt.title('Interaction matrix')
        plt.colorbar(plt.imshow(mat, cmap='seismic_r', vmin=-np.max(np.abs(mat)), vmax=np.max(np.abs(mat))
                 ),ax=plt.gca())
        plt.subplot(133)
        offdiag=(np.eye(mat.shape[0])==0)&(np.isnan(corr)==0)
        xs,ys=mat[offdiag],corr[offdiag]
        if xs.any():
            strength=np.mean(ys)
            plt.scatter(xs,ys)
            pxs, pys, desc, slope, intercept, r, p, stderr = linfit(xs,ys)
            plt.plot(pxs,pys,lw=2,color='k')
            plt.xlabel('Interaction matrix'),plt.ylabel('Spatial correlation')
            plt.title(r'{:.2f} + {:.2f} x ($R^2$={:.2f}, p={:.1g})'.format(intercept,slope,r**2,p))
        else:
            r,slope,strength=0,0,0
        dic.update({'checker_r2':r**2,'checker_slope':slope ,'checker_strength':strength})

        if save:
            for f in figs:
                figs[f].savefig(fpath+'{}.png'.format(f) )
            plt.close('all')
        else:
            plt.show()
    # code_debugger()
    return table


def summary_plots(path,axes=None,save=0,values=None,**kwargs):
    """Creates summary measurements for each simulation, store them all in summary.csv, along with all simulation parameters."""

    df=None
    if not kwargs.get('rerun',0):
        try:
            df=pd.read_json(path+'summary.csv')
        except:
            pass
    if df is None:
        if kwargs.get('detailed',1):
            table=detailed_plots(path,save,**kwargs)
        else:
            # Measure summary quantities
            files=pd.read_csv(path+'files.csv')
            table=[]
            for idx,f in files.iterrows():
                model=LandscapeModel.load(f['path'])
                dic = dict(model.export_params())
                dic['path']=f['path']
                basic_measures(model,dic)
                table.append(dic)
        df=pd.DataFrame(table)
        df.to_json(path+'summary.csv')

    if not axes is None:
        axes = [a if not hasattr(a, '__iter__') else a[0] for a in axes]
        axes=[a for a in axes if  a !='sys']
    else:
        model = LandscapeModel.load(df['path'].values[0])
        refprm=model.export_params()
        axes=[k for k in refprm  if len(np.unique(df[k].values))>1  and k != 'sys']

    figs={}
    try:
        import seaborn as sns
        def dist_axes(axes,names=None,max=3):
            if not names:
                names=['col','row','hue','style']
            names=names[:max]
            return dict(zip(names[:len(axes)],axes[:len(names)]) )
        if 'SAR' in df:
            sardf=pd.DataFrame([dict(zip(axes,row[axes].values)+[('area',x),('SAR',y)] ) for i,row in  df.iterrows() for x,y in zip(*row['SAR'])])
            # code_debugger()
            figs['SAR']=sns.relplot(x='area',y='SAR',data=sardf,kind="line",**dist_axes(axes,max=4))
            plt.xscale('log')


        # Figure: effects of dispersal
        results=[ ('alpha_div', 'Alpha diversity'), ('alpha_div_min', 'Alpha diversity min'),
                  ('alpha_div_max', 'Alpha diversity max'),('biomass_tot','Total biomass'),
                  ] + list(values)
        for ax in axes:
            df[ax]= [np.round(z,int(max(1,np.ceil(-np.log10(z)))))  for z in df[ax].values ]
        for value,title in results:
            # panel, ax = auto_subplot(panel, len(results))
            axd=[ax for ax in axes if not 'dispersal' in ax]
            figs[value]=sns.catplot(x='dispersal_mean', y=value,  data=df,kind='box',
                                         height=6, aspect=.75, linewidth=2.5, **dist_axes(axd))
            plt.suptitle(title)
    except Exception as e:
        print e

    if save:
        if isinstance(save,basestring) or isinstance(save,Path):
            outpath=save
        else:
            outpath = path + Path('plots')
        for f in figs:
            figs[f].savefig(outpath+'{}.png'.format(f) )
    else:
        plt.show()
