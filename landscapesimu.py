from landscapemodel import *
from plots import *


# ============= SIMULATIONS =================

def loop(axes,path,rerun=0,multiscale=1,tmax=50,tsample=10,keep='all',**kwargs):
    """Iterate over parameters given in axes.
    Examples:
        axes=[ ('species', [1,10,100] ), ('landx',[10,100])   ]

        will iterate over all possible combinations of species number and landscape width among given values,
        save results in directory tree in path, and plot results.
    """

    if not 'files.csv' in os.listdir(path):
        rerun=1

    # Get default parameter values (see landscapemodel.py for parameter list)
    prm=deepcopy(LandscapeModel.dft_prm)

    for i in prm:
        if 'multiscale' in str(prm[i]):
            if multiscale=='all' or i in multiscale:
                prm[i]['multiscale']=1

    # Update parameters
    for i,j in kwargs.items():
        if i in prm:
            prm[i]=j

    if rerun:
        Looper(LandscapeModel).loop(tmax=tmax,tsample=tsample,keep=keep,path=path,axes=axes,parameters=prm)
        rebuild_filelist(path)


# Only run the following code if this script is executed directly
if __name__=='__main__':

    # Folder containing simulations
    path = Path('RESULTS/beflandscape')

    # Multiscale?
    multiscale=''#'all'

    # Simulation length
    tmax=500

    # Get options from command line
    for i in sys.argv:
        if 'path=' in i:
            path=Path('RESULTS/{}'.format(i.split('=')[-1]) )
        if 'multiscale' in i:
            if '=' in i:
                multiscale=i.split('=')[-1]
            else:
                multiscale='all'
        if 'tmax=' in i :
            tmax=ast.literal_eval(i.split('=')[-1])

    # Loop over parameters given in axes
    # 'sys' is a dummy parameter allowing for replicas
    axes = [('sys', [0])]
    loop(axes=axes,path=path,tmax=tmax,rerun='rerun' in sys.argv,multiscale=multiscale)

    # Plots
    # Detailed plots for each simulation
    detailed_plots(path)

    # Summary plots for BEF over parameter values
    summary_plots(path,axes=axes)