from landscapemodel import *
from plots import *

def FTtest(path='TEST/FTtest',species=6,tmax=1):
    """Test Fourier Transform optimization"""
    path=Path(path).mkdir()
    if not 'files.csv' in os.listdir(path):
        rerun=1

    prm=deepcopy(LandscapeModel.dft_prm)
    prm['species']=species
    prm['landx']=prm['landy']=32
    prm['dispersal']['mean']=1

    model=LandscapeModel(parameters=prm)
    model.evol(use_Fourier=0,dt=0.01,tmax=tmax)
    direct =deepcopy(model.results)
    model.evol(use_Fourier=1,reseed=0,init='restart',dt=0.01,tmax=tmax)
    fourier =deepcopy(model.results)


    for s in range(species):
        plt.figure()
        plt.suptitle('Species {}'.format(s) )
        panel=0
        for i in [0,-1]:#range(len(direct['n']))
            panel, ax = auto_subplot(panel, 4)
            d,f=direct['n'][i][s],fourier['n'][i][s]
            span=np.min(np.minimum(d,f)),np.max(np.maximum(d,f))
            plt.colorbar(plt.imshow(d,vmin=span[0],vmax=span[1]),ax=ax)
            t=direct['t'][i]
            plt.title('Direct t={}'.format(t) )
            panel, ax = auto_subplot(panel, 4)
            plt.colorbar(plt.imshow(f,vmin=span[0],vmax=span[1]),ax=ax)
            plt.title('Fourier t={}'.format(t) )

    plt.show()
    code_debugger()

if __name__=='__main__':
    FTtest()