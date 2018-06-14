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

def reorder_FT(x,k=None):
    """Fourier transforms are typically in order [0..k(N-1),-k(N)...-k(1)].
    This reorders them to [-k(N)..k(N-1)] or back."""
    if k is None:
        lx, ly = x.shape[-2:]
        kx,ky=fft.fftfreq(lx), fft.fftfreq(ly)
    else:
        kx,ky=k
    x2=x[np.ix_(*[range(z) for z in x.shape[:-2]]+[ np.argsort(kx),np.argsort(ky) ])]
    return x2

def checkpos_FT(x,idx=0):
    """Checks positivity of first element of x, by looking at conjugate symmetry of Fourier Transform."""
    t=x
    while len(t.shape)>2:
        t=t[idx]
    lx,ly=t.shape
    check=np.max(np.abs(t[lx / 2:,1: ] - np.conj(t[1:lx / 2 + 1,1:][::-1, ::-1])))
     #np.max(np.abs(t[lx / 2:, ly / 2:] - np.conj(t[1:lx / 2 + 1, 1:ly / 2 + 1][::-1, ::-1])))
    # code_debugger()
    return check

def generate_noise(shape,method='fft', **dprm):
    """Generate noisy landscape (for now only method=fft and filter works convincingly"""
    samples = dprm.get('samples', 500)
    dimres = []
    if method=='fft':
        ft = np.random.random(shape)
        ft-=np.mean(ft)
        ft[0,0]=10
        if len(shape)==2:
            color=dprm.get('color',1)
            lx, ly = ft.shape
            kx, ky = fft.fftfreq(lx), fft.fftfreq(ly)
            k=np.add.outer(kx ** 2, ky ** 2)**0.5
            ft=ft / (0.00000001+k** (color )) * np.exp(-k/np.max(k)*1./dprm.get('cutoff',1000))
            res=fft.fft2(ft )
            res=reorder_FT(res.real+res.imag)
            plt.imshow(res)
            plt.show()
        else:
            raise Exception('Not implemented yet')
    elif method=='filter':
        res = np.random.random(shape)
        res = ndimage.gaussian_filter(res, sigma=5)
    elif method=='direct':
        print 'Experimental, failed'
        for dim in range(len(shape)):
            freqs = np.logspace(0, np.log(shape[dim] / 100.), samples)
            amps = freqs ** dprm.get('color', 1)
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
