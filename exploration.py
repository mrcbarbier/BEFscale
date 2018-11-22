from landscapesimu import *

# model=LandscapeModel(parameters=prm)

def single_species(path,tmax=1000,nsample=50,rerun=0,replot=0,**kwargs):
    prm=deepcopy(LandscapeModel.dft_prm)
    prm['species']=1
    prm['death']=10**-2
    prm['landx']=prm['landy']=64
    prm['dispersal']['mean']=1
    prm['n'].update({'range':(1.,2.), 'spatial':'spotty',
                     'distribution_prm':{'axes':(1,2), 'number':1,'shape':'blot', 'width':.2 }})
    prm['envniche']['width']['range']=(5,6)

    axes=[
        ('environment_cutoff', [0.05, .2]),  # Inverse of environment spatial scale
        ('dispersal_mean', [0.01, .2]),#,10.]),  #
        ('sys', range(1)),  # Dummy variable
        ('n_spatial', ['spotty',]),# 'uniform', ]),  # Mean dispersal strength
        ('competition_scale', [0.01, 5.]),  #
    ]


    loop(axes=axes,path=path,tmax=tmax,nsample=nsample,samplescale='linear',rerun=rerun,
         reseed=0,use_Fourier=0,method='scipy',dft_prm=prm,**kwargs)

    summary_plots(path,axes=axes,save=1,movie=0,detailed=1,rerun=rerun or replot)


if __name__=='__main__':
    path = Path('RESULTS/single_species')     # Folder containing simulations
    single_species(path,rerun='rerun' in sys.argv ,replot= 'replot' in sys.argv)