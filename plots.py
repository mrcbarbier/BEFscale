from landscapemodel import *
import matplotlib.figure as mpfig


# ============= PLOTS =================

def detailed_plots(path,save=0):
    """Plots with detailed information on each simulation outcome"""


    files=pd.read_csv(path+'files.csv')
    for idx,f in files.iterrows():
        model=LandscapeModel.load(f['path'])
        prm,data=model.prm,model.data
        death=prm['death']

        Nf=model.results['n'][-1]/death
        N=prm['species']

        #Plots of initial species abundance per patch
        figNini=plt.figure(figsize=np.array(mpfig.rcParams['figure.figsize'])*(3,3) )
        panel=0
        for i in range(N):
            panel,ax=auto_subplot(panel,N)
            plt.imshow(model.results['n'][0][i])#,vmin=0,vmax=1.2)
            plt.colorbar()
            plt.title('Species {}'.format(i))
        plt.suptitle('Initial abundance')

        #Plots of species abundance per patch
        figNfin=plt.figure(figsize=np.array(mpfig.rcParams['figure.figsize'])*(3,3)  )
        panel=0
        for i in range(N):
            panel,ax=auto_subplot(panel,N)
            plt.imshow(Nf[i])#,vmin=0,vmax=1.2)
            plt.colorbar()
            plt.title('Species {}'.format(i))
        plt.suptitle('Final abundance')

        #Species environmental niches
        figenv=plt.figure()
        panel=0
        panel, ax = auto_subplot(panel, 2)
        plt.imshow(data['environment'])
        plt.title('Environment value')
        plt.colorbar()
        panel, ax = auto_subplot(panel, 2)
        plt.errorbar(np.arange(N),data['envniche_pos'],yerr=data['envniche_width'] )
        plt.title('Species abiotic niche')

        #Role of different interactions
        figdynamics=plt.figure(figsize=np.array(mpfig.rcParams['figure.figsize'])*(3,3) )
        panel=0
        dx,typfluxes,fluxes=model.get_dx(None, Nf, calc_fluxes=1)
        for i in range(N):
            panel,ax=auto_subplot(panel,N)
            for t, typ in enumerate(typfluxes):
                ff=fluxes[t,i]/np.abs(dx[i])
                ff=np.mean(ff)  #Average over all patches
                # f=f.ravel()[np.argmax(Nf[i])]  #Values at most abundant patch
                plt.bar(t,ff)
            plt.title('Species {}'.format(i))
        plt.legend(typfluxes)
        plt.suptitle('Dynamical effect of interactions')


        #Biomass pyramid
        plt.figure()
        plt.scatter(data['size'],np.sum(Nf,axis=(1,2)))
        plt.xlabel('Size')
        plt.yscale('log')
        plt.ylabel('Biomass')

        #BEF
        figBEF=plt.figure()
        B=np.sum(Nf,axis=0)  #Total Biomass
        D=np.sum(Nf>death,axis=0) #Diversity
        n=Nf/B
        DS=-np.sum(n*np.log(np.clip(n,10**-15,None) ),axis=0)  #Shannon diversity

        patches=pd.DataFrame([dict( zip(('B','n','D','DShannon'), x)  ) for x in zip(B.ravel(),n.ravel(),D.ravel(),DS.ravel()) ]  )
        plt.subplot(121)

        # Group by species diversity
        patchesD=patches.groupby('D').median()
        plt.scatter(D.ravel(),B.ravel())
        #plt.yscale('log')
        #plt.xscale('log')
        plt.plot(patchesD.index,patchesD['B'],lw=2,c='k' )
        plt.title('Diversity')
        plt.subplot(122)

        # Group by Shannon
        bins = np.linspace(patches['DShannon'].min(), patches['DShannon'].max(), 5)
        patchesDS=patches.groupby(np.digitize(patches['DShannon'], bins)).median()
        plt.scatter(DS.ravel(),B.ravel())
        # plt.yscale('log')
        # plt.xscale('log')
        plt.plot(patchesDS['DShannon'],patchesDS['B'],lw=2,c='k' )
        plt.title('Shannon diversity')

        if save:
            if isinstance(save,basestring) or isinstance(save,Path):
                outpath=save
            else:
                outpath = path + Path('plots')
            fpath=Path(outpath+'_'.join([n for n in Path(f['path']).split() if not n in Path(outpath).split() ]))
            fpath.mkdir()
            figNini.savefig(fpath+'N_initial.png')
            figNfin.savefig(fpath+'N_final.png')
            figenv.savefig(fpath+'Environment.png')
            figdynamics.savefig(fpath+'dynamics.png')
            figBEF.savefig(fpath+'BEF.png')
        else:
            plt.show()

        # code_debugger()


def measures(model,dic):
    """Measure various summary quantities for a given model, store them in dictionary dic."""
    Nf = model.results['n'][-1]
    B=np.sum(Nf,axis=0)  #Total Biomass per patch
    n=Nf/B
    B = np.sum(Nf)  # Total Biomass on landscape
    D = np.sum(Nf > 0.01, axis=0)  # Diversity
    DS = -np.sum(n * np.log(np.clip(n, 10 ** -15, None)), axis=0)  # Shannon diversity
    dic.update({ 'biomass_tot':np.sum(B),'alpha_div':np.mean(D),'alpha_shannon':np.mean(DS)})


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
