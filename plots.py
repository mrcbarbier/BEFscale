from landscapemodel import *
import matplotlib.figure as mpfig




def measures(model,dic):
    """Measure various quantities for a given model, store summary in dictionary dic."""
    Nf = model.results['n'][-1]
    growth=model.data['growth']

    B=np.sum(Nf,axis=0)  #Total Biomass per patch
    basal=np.where(np.max(growth,axis=(1,2))>0 )[0] #Basal species
    Bbasal=np.sum(Nf[basal ],axis=0)  #Total producer biomass per patch
    prod=np.sum(growth*Nf,axis=0)  #Total production per patch

    alive=(Nf > model.prm['death'])
    n=Nf/B # Relative biomass
    D = np.sum(alive, axis=0)  # Diversity
    DS = -np.sum(n * np.log(np.clip(n, 10 ** -15, None)), axis=0)  # Shannon diversity
    dic.update({ 'biomass_tot':np.sum(B),'biomass_basal':np.sum(Bbasal),
                 'production_tot':np.sum(prod),
                 'alpha_div':np.mean(D),'alpha_shannon':np.mean(DS),
                 'gamma_div':np.sum(np.sum(alive,axis=(1,2))>0 ) })

    return B, Bbasal, D, DS

# ============= PLOTS =================

def detailed_plots(path,save=0,movie=0):
    """Plots with detailed information on each simulation outcome"""

    if save:
        if isinstance(save, basestring) or isinstance(save, Path):
            outpath = save
        else:
            outpath = path + Path('plots')

    files=pd.read_csv(path+'files.csv')
    for idx,f in files.iterrows():
        model=LandscapeModel.load(f['path'])
        prm,data=model.prm,model.data
        death=prm['death']

        if save:
            fpath=Path(outpath+'_'.join([n for n in Path(f['path']).split() if not n in Path(outpath).split() ]))
            fpath.mkdir()

        nf=model.results['n'][-1]
        N=prm['species']

        #Species trajectories
        figtraj=plt.figure()
        nsum=np.sum(model.results['n']/death, axis=(2, 3))
        for i in zip(*nsum):
            plt.plot(model.results['t'],i)
        plt.yscale('log')
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
                figNini=fig
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

        figNfin=fig
        plt.suptitle('Final abundance')


        #Species environmental niches
        figenv=plt.figure()
        panel=0
        panel, ax = auto_subplot(panel, 2)
        plt.imshow(data['environment'])
        plt.title('Environment value')
        plt.colorbar()
        panel, ax = auto_subplot(panel, 2)

        c,w=data['envniche_pos'], data['envniche_width']
        plt.vlines(np.arange(N),c-w,c+w )
        plt.title('Species abiotic niche')


        #Food web
        if model.prm['trophic']['ON']:
            figweb=plt.figure(figsize=np.array(mpfig.rcParams['figure.figsize'])*(2,2) )
            draw_network(model.data['trophic'],ypos=dict(zip(range(N),model.data['size'] )),newfig=0,hold=1)
            plt.ylabel('Body size')
            plt.title('Food web')
        else:
            figweb=None

        #Role of different interactions
        figdynamics=plt.figure(figsize=np.array(mpfig.rcParams['figure.figsize'])*(3,3) )
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
        figpyr=plt.figure()
        plt.scatter(data['size'],np.sum(nf/death,axis=(1,2)))
        plt.xlabel('Size')
        plt.yscale('log')
        plt.ylabel('Abundance')
        plt.title('Size spectrum')

        #SAR
        figSAR=plt.figure()
        nstart=100
        startpos=(np.random.randint(prm['landx'],size=nstart),np.random.randint(prm['landy'],size=nstart) )
        dist=(np.add.outer(startpos[0],-startpos[0])**2+np.add.outer(startpos[1],-startpos[1])**2)
        SARS=[[],[]]
        for i in range(nstart):
            order=np.argsort(dist[i])
            xs=dist[i,order]
            ys=[nf[:,startpos[0][o],startpos[1][o]]>death for o in order]
            ys=np.cumsum(ys,axis=0)
            ys=np.sum(ys>0,axis=1)
            plt.plot(xs,ys)
            SARS[0]+=list(xs)
            SARS[1]+=list(ys)
            # code_debugger()
        xs,ys=np.array(SARS)
        bins=np.logspace(np.min(np.log10(xs[xs>0])),np.max(np.log10(xs[xs>0])),20)
        ibins=np.digitize(xs,bins)
        plt.plot(bins,[np.mean(ys[ibins==i] )  for i in range(len(bins)) ] ,lw=3,c='k')
        plt.xscale('log')
        plt.title('SAR')
        plt.xlabel('Area')
        plt.ylabel('Diversity')


        #BEF
        figBEF=plt.figure()
        summaries={}
        B,Bbasal,D,DS=measures(model,summaries)

        patches=pd.DataFrame([dict( zip(('B','Bbasal','D','DShannon'), x)  ) for x in zip(B.ravel(),Bbasal.ravel(),
                                                                                          D.ravel(),DS.ravel()) ]  )

        panel=0
        for name,val in [('Biomass','B'),('Basal biomass','Bbasal')]:
            panel,ax=auto_subplot(panel,4)
            # Group by species diversity
            patchesD=patches.groupby('D').median()
            plt.scatter(patches['D'],patches[val])
            #plt.yscale('log'),plt.xscale('log')
            plt.plot(patchesD.index,patchesD[val],lw=2,c='k' )
            plt.xlabel('Diversity')
            plt.ylabel(name)
            panel,ax=auto_subplot(panel,4)

            # Group by Shannon
            bins = np.linspace(patches['DShannon'].min(), patches['DShannon'].max(), 5)
            patchesDS=patches.groupby(np.digitize(patches['DShannon'], bins)).median()
            plt.scatter(patches['DShannon'],patches[val])
            # plt.yscale('log'), plt.xscale('log')
            plt.plot(patchesDS['DShannon'],patchesDS[val],lw=2,c='k' )
            plt.xlabel('Shannon diversity')
            plt.ylabel(name)

        if save:
            figNini.savefig(fpath+'N_initial.png')
            figNfin.savefig(fpath+'N_final.png')
            figenv.savefig(fpath+'Environment.png')
            figdynamics.savefig(fpath+'dynamics.png')
            figBEF.savefig(fpath+'BEF.png')
            figSAR.savefig(fpath+'SAR.png')
            figtraj.savefig(fpath+'N_traj.png')
            figpyr.savefig(fpath+'N_pyramid.png')
            if not figweb is None:
                figweb.savefig(fpath+'foodweb.png')
        else:
            plt.show()

        # code_debugger()



def summary_plots(path,axes=None,save=0):
    """Creates summary measurements for each simulation, store them all in summary.csv, along with all simulation parameters."""

    # Measure summary quantities
    files=pd.read_csv(path+'files.csv')
    table=[]
    for idx,f in files.iterrows():
        model=LandscapeModel.load(f['path'])
        dic={'path':f['path'] }
        dic.update( model.export_params())
        measures(model,dic)
        table.append(dic)
    df=pd.DataFrame(table)
    df.to_csv('summary.csv')

    # Figure: effects of dispersal
    figdisp=plt.figure()
    panel = 0
    results=[ ('alpha_div', 'Alpha diversity'), ('biomass_tot','Total biomass')  ]
    for value,title in results:
        panel, ax = auto_subplot(panel, len(results))
        plt.plot(df['dispersal_mean'],df[value])
        plt.xlabel('Dispersal')
        plt.title(title)

    if save:
        if isinstance(save,basestring) or isinstance(save,Path):
            outpath=save
        else:
            outpath = path + Path('plots')
        figdisp.savefig(outpath+'dispersal.png')
    else:
        plt.show()
