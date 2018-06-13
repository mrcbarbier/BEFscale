import numpy as np, pandas as pd,scipy.integrate as scint, scipy.ndimage as ndimage, itertools
import time,sys, matplotlib.pyplot as plt, scipy.fftpack as fft,scipy.signal as ssignal, os, ast
from copy import deepcopy


def code_debugger(skip=0):
    """Utility: """
    import code
    import inspect
    stack=inspect.stack()
    dic = {}
    dic.update(stack[1+skip][0].f_globals)
    dic.update(stack[1+skip][0].f_locals)
    code.interact(local=dic)

def generate_noise(shape,method='filter', **dprm):
    """Generate noisy landscape (experimental, for now only method=filter works convincingly"""
    samples = dprm.get('samples', 500)
    dimres = []
    if method=='filter':
        res = np.random.random(shape)
        res = ndimage.gaussian_filter(res, sigma=5)
        # res+=ndimage.gaussian_filter(res, sigma=80)#*np.max(res)
    elif method=='direct':
        for dim in range(len(shape)):
            freqs = np.logspace(0, np.log(shape[dim] / 100.), samples)
            amps = freqs ** dprm.get('spectralexp', 0)
            phase = np.random.random(samples)
            xs = np.zeros(shape)
            dx = np.linspace(0, 1, shape[dim]).reshape(
                [1 for z in range(0, dim)] + [shape[dim]] + [1 for z in range(dim + 1, len(shape))])
            ps = np.multiply.outer(*[np.linspace(0, 1, sh) for sh in shape])
            xs = xs + dx
            dimres.append((np.sum([a * np.exp(1j * 2 * np.pi * (xs * f + p * np.cos(2 * np.pi * ps))) for f, a, p in
                                   zip(freqs, amps, phase)], axis=0)))
        res = np.real(np.multiply(*dimres))

    else:
        freqs = np.logspace(0, 2, samples)
        amps = freqs ** dprm.get('spectralexp', 0)
        phase = np.random.random((samples, len(shape)))
        xs = [np.zeros(shape) + np.linspace(0, 1, sh).reshape(
            [1 for z in range(0, dim)] + [sh] + [1 for z in range(dim + 1, len(shape))]) for dim, sh in
              enumerate(shape)]
        res = np.real(np.sum([a * np.exp(1j * 2 * np.pi * np.add(*[x * f + pp for x, pp in zip(xs, p)])) for f, a, p in
                              zip(freqs, amps, phase)], axis=0))
    return res

def dumps(fil,obj):
    """Write object to file"""
    fil.write(str(obj))

def loads(fil):
    """Load object from file"""
    txt=''
    for l in fil:
        txt+=l.strip()
    return eval(txt,{},{'array':np.array,'nan':np.nan})


def rebuild_filelist(path,verbose=True):
    """Run through directory tree to look for model results to add to files.csv"""
    if verbose:
        print 'Rebuilding files.csv for {}'.format(path)
    final=None
    idx=0
    for root,dirs,files in os.walk(str(path)):
        if str(root)!=str(path):
            if 'files.csv' in files:
                pass
            elif 'model.dat' in files and 'results.csv' in files:
                print 'Missing files.csv, regenerating'
                dic={'path':root }
                #dic.update( eval(open('model.dat','r' ) ) )
                df=pd.DataFrame([dic])
                df.to_csv(Path(root)+'files.csv')
            else:
                continue
            if verbose:
                print '{} Found files.csv in {}'.format(idx, root)
                idx+=1
            if final is None:
                final=pd.read_csv(Path(root)+'files.csv',index_col=0)
            else:
                final=final.append(pd.read_csv(Path(root)+'files.csv',index_col=0),ignore_index=1)
            final.loc[final.index[-1],'path']=str(Path(root).osnorm())

    if not final is None:
        final.to_csv(path+'files.csv')
    else:
        print 'Error: No file found! {}'.format(path)


class Path(str):
    '''Strings that represent filesystem paths.
    Overloads __add__:
     - when paths are added, gives a path
     - when a string is added, gives a string'''
    def __add__(self,x):
        import os
        if isinstance(x,Path):
            return Path(os.path.normpath(os.path.join(str(self),x)))
        return os.path.normpath(os.path.join(str(self),x))

    def norm(self):
        import os
        return Path(os.path.normpath(str(self)))

    def osnorm(self):
        """Deal with different separators between OSes."""
        import os
        if os.sep=='/' and "\\" in str(self):
            return Path(os.path.normpath(str(self).replace('\\','/' )))
        elif os.sep=='\\' and "/" in str(self):
            return Path(os.path.normpath(str(self).replace('/','\\' )))
        else:
            return self.norm()

    def prev(self):
        import os
        lst=self.split()
        path=os.path.join(lst[:-1])
        return path.osnorm()

    def split(self,*args):
        """"""
        import os
        lst=[]
        cur=os.path.split(self.norm())
        while cur[-1]!='':
            lst.insert(0,cur[-1])
            cur=os.path.split(cur[0])
        return lst

    def mkdir(self,rmdir=False):
        """Make directories in path that don't exist. If rmdir, first clean up."""
        import os
        if rmdir:
            os.rmdir(str(self))
        cur=Path('./')
        for intdir in self.split():
            cur+=Path(intdir)
            if not os.path.isdir(cur):
                os.mkdir(cur)
        return self

    def copy(self):
        return Path(self)

    def strip(self,*args):
        '''Return string without final / or \\ to suffix/modify it.'''
        return str(self).strip('\/')


def auto_subplot(panel,nbpanels,rows=None,projection=None,return_all=0):
    """Helps deal with pyplot's subplots, automatically incrementing panel index."""
    i=panel
    if rows is None:
        panels_per_row=np.ceil(np.sqrt(nbpanels) )
    else:
        panels_per_row=np.ceil(nbpanels/rows).astype('int')
    nbrows=np.ceil(nbpanels*1./panels_per_row)
    ax=plt.subplot(nbrows,panels_per_row,i +1 ,projection=projection)
    panel+=1
    if return_all:
        return ax,nbrows,panels_per_row
    return panel,ax
