from landscapemodel import *


# ============= PLOTS =================

def detailed_plots(path):
    """Plots with detailed information on each simulation outcome"""

    files=pd.read_csv(path+'files.csv')
    for idx,f in files.iterrows():
        model=LandscapeModel.load(f['path'])
        Nf=model.results['n'][-1]
        prm,data=model.prm,model.data
        N=prm['species']
        death=prm.get('death',10**-15) #Death threshold


        #Plots of species abundance per patch
        plt.figure()
        panel=0
        for i in range(N):
            panel,ax=auto_subplot(panel,N)
            plt.imshow(Nf[i])#,vmin=0,vmax=1.2)
            plt.colorbar()
            plt.title('Species {}'.format(i))
        plt.suptitle('Abundance')

        #Species environmental mortality
        plt.figure()
        panel=0
        for i in range(N):
            panel,ax=auto_subplot(panel,N)
            plt.imshow(data['mortality'][i])
            plt.title('Species {}'.format(i))
        plt.suptitle('Mortality')

        #Role of different interactions
        plt.figure()
        panel=0
        dx,typfluxes,fluxes=model.get_dx(None, Nf, calc_fluxes=1)
        for i in range(N):
            panel,ax=auto_subplot(panel,N)
            for t, typ in enumerate(typfluxes):
                f=fluxes[t,i]/np.abs(dx[i])

                f=np.mean(f)  #Average over all patches
                # f=f.ravel()[np.argmax(Nf[i])]  #Values at most abundant patch
                plt.bar(t,f)
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
        plt.figure()
        B=np.sum(Nf,axis=0)  #Total Biomass
        D=np.sum(Nf>0.01,axis=0) #Diversity
        n=Nf/B
        D2=-np.sum(n*np.log(np.clip(n,10**-15,None) ),axis=0)  #Shannon diversity

        patches=pd.DataFrame([dict( zip(('B','n','D','DShannon'), x)  ) for x in zip(B.ravel(),n.ravel(),D.ravel(),D2.ravel()) ]  )
        plt.subplot(121)
        patchesD=patches.groupby('D').median()
        plt.scatter(D.ravel(),B.ravel())
        plt.yscale('log')
        plt.xscale('log')
        plt.plot(patchesD.index,patchesD['B'],lw=2,c='k' )
        plt.title('Diversity')
        plt.subplot(122)
        bins = np.linspace(patches['DShannon'].min(), patches['DShannon'].max(), 5)
        patchesDS=patches.groupby(np.digitize(patches['DShannon'], bins)).median()
        plt.scatter(D2.ravel(),B.ravel())
        plt.yscale('log')
        plt.xscale('log')
        plt.plot(patchesDS['DShannon'],patchesDS['B'],lw=2,c='k' )
        plt.title('Shannon diversity')
        plt.show()

        # code_debugger()



# ============= SIMULATIONS =================

def loop(path,rerun=0,**kwargs):
    prm=deepcopy(LandscapeModel.dft_prm)
    for i,j in kwargs.items():
        if i in prm:
            prm[i]=j
    axes=[('sys',[0]) ]
    if rerun:
        Looper(LandscapeModel).loop(tmax=500,tsample=10,keep='all',path=path,axes=axes,parameters=prm)
        rebuild_filelist(path)

if __name__=='__main__':

    # Folder containing simulations
    path = Path('beflandscape')

    loop(path=path,rerun='rerun' in sys.argv)

    # Plots
    detailed_plots(path)