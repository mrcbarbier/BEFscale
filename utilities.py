import numpy as np, pandas as pd,scipy.integrate as scint, scipy.ndimage as ndimage, itertools
import time,sys, matplotlib.pyplot as plt, scipy.fftpack as fft,scipy.signal as ssignal, os
from copy import deepcopy


def code_debugger(skip=0):
    import code
    import inspect
    stack=inspect.stack()
    dic = {}
    dic.update(stack[1+skip][0].f_globals)
    dic.update(stack[1+skip][0].f_locals)
    code.interact(local=dic)

def dumps(fil,obj):
    fil.write(str(obj))

def loads(fil):
    txt=''
    for l in fil:
        txt+=l.strip()
    return eval(txt,{},{'array':np.array,'nan':np.nan})


def rebuild_filelist(path,verbose=True):
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
            final.set_value(final.index[-1],'path',str(Path(root).osnorm()))

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

    def copy(self):
        return Path(self)

    def strip(self,*args):
        '''Return string without final / or \\ to suffix/modify it.'''
        return str(self).strip('\/')


def auto_subplot(panel,nbpanels,rows=None,projection=None,return_all=0):
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
