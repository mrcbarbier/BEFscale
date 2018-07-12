from landscapemodel import *
from plots import *


def ABtest(comparison,path,species=1,tmax=100,tshow=(0,-1),**kwargs):
    """A/B testing: plot results side-by-side for sets of options in comparison"""
    path=Path(path).mkdir()
    prm=deepcopy(LandscapeModel.dft_prm)
    prm['species']=species
    prm['landx']=prm['landy']=32
    prm['dispersal']['mean']=1

    model=LandscapeModel(parameters=prm)
    results=[]
    for comp in comparison:
        dic={}
        dic.update(kwargs)
        if results:
            dic['reseed']=0
            dic['init']='restart'
        dic.update(comp)
        model.evol(dt=0.01,tmax=tmax,**dic)
        results.append(deepcopy(model.results))

    for s in range(species):
        plt.figure()
        plt.suptitle('Species {}'.format(s) )
        panel=0
        for i in tshow:
            span=min([np.min(r['n'][i][s]) for r in results ]) , max([np.max(r['n'][i][s]) for r in results ])
            for idx,comp, res in zip(range(len(results)),comparison,results):
                panel, ax = auto_subplot(panel, len(tshow) * len(results))
                plt.colorbar(plt.imshow(res['n'][i][s],vmin=span[0],vmax=span[1]),ax=ax)
                t=res['t'][i]
                plt.title('{} t={}'.format(comp.get('title',idx), t) )

    if kwargs.get('debug',0):
        code_debugger()
    else:
        plt.show()

def FTtest(path='TEST/FTtest',**kwargs):
    """Test Fourier Transform optimization"""
    comparison=[{'use_Fourier':0,'title':'Direct'},{'use_Fourier':1,'title':'Fourier'}  ]
    ABtest(comparison,path,**kwargs)

def algotest(path='TEST/algotest',**kwargs):
    """Test dependence on integration algorithm"""
    comparison=[{'method':'Euler','title':'Euler'},{'method':'scipy','title':'scipy+Fourier','use_Fourier':1}  ]
    ABtest(comparison,path,**kwargs)


if __name__=='__main__':
    # FTtest()
    algotest(debug='debug' in sys.argv)