from landscapemodel import *



# ============= SIMULATIONS =================

def loop(axes,path,rerun=0,multiscale=1,**kwargs):
    """Iterate over parameters given in axes.
    Examples:
        axes=[ ('species', [1,10,100] ), ('landx',[10,100])   ]

        will iterate over all possible combinations of species number and landscape width among given values.
        """



    # Get default parameter values (see landscapemodel.py for parameter list)
    prm=deepcopy(Landscape.dft_prm)

    for i in prm:
        if 'multiscale' in prm[i]:
            if multiscale=='all' or i in multiscale:
                prm[i]['multiscale']=1

    # Update parameters
    for i,j in kwargs.items():
        if i in prm:
            prm[i]=j

    if rerun:
        Looper(LandscapeModel).loop(tmax=500,tsample=10,keep='all',path=path,axes=axes,parameters=prm)
        rebuild_filelist(path)


# Only run the following code if this script is executed directly
if __name__=='__main__':

    # Folder containing simulations
    path = Path('RESULTS/beflandscape')

    # Multiscale?
    multiscale='all'

    # Get options from command line
    for i in sys.argv:
        if 'path=' in i:
            path=Path('RESULTS/{}'.format(i.split('=')[-1]) )
        if 'multiscale' in i:
            if '=' in i:
                multiscale=i.split('=')[-1]
            else:
                multiscale='all'

    # Loop over parameters given in axes
    # 'sys' is a dummy parameter allowing for replicas
    axes = [('sys', [0])]
    loop(axes=axes,path=path,rerun='rerun' in sys.argv,multiscale=multiscale)

    # Plots
    detailed_plots(path)