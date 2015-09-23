

__all__ = []

import gzip

from tpr   import TPR
from trr   import TRR
from conf  import CONF,PDB,GRO
from top   import TOP
from index import Index


class GMXIOError(Exception):
    def __init__(self,msg): self.msg = msg
    def __str__(self):      return self.msg


_loaders = [
    (".pdb",    PDB),
    (".pdb.gz", PDB),
    (".ent",    PDB),
    (".ent.gz", PDB),
    (".gro",    GRO),
    (".gro.gz", GRO),
    (".trr",    TRR),
    (".tpr",    TPR),
    (".ndx",    Index),
    (".ndx.gz", Index),
    (".select", Index),
    (".top",    TOP),
    (".top.gz", TOP)
    ]


_what = {
    "ndx": Index,
    "trr": TRR,
    "gro": GRO,
    "pdb": PDB,
    "tpr": TPR,
    "top": TOP
    }
   
 
def open(stream,**kwargs):
    """Open the file using the loader determined from the extension"""

    if not stream:
        return None

    what = kwargs.get("what")
    fun  = None

    if type(stream) == str:
        if what:
            fun = _what.get(what)
            if not fun:
                raise GMXIOError("Unknown file type for %s: %s"%(stream,what))
        else:
            for ext,fn in _loaders:
                if stream.endswith(ext):
                    fun = fn
            if not fun:
                raise GMXIOError("Unable to determine file type from extension for %s"%stream)
    
    if not fun and not what:
        # Only PDB/GRO can be handled like this 
        fun = CONF
    
    if not fun:
        raise GMXIOError("Unknown type for open stream: %s"%what)

    # Extract keywords relevant to function
    
    if hasattr(fun,"func_code"):
        # Function arguments
        args = fun.func_code.co_varnames
    else:
        # Constructor arguments, excluding 'self'
        args = fun.__init__.func_code.co_varnames[1:]
    rubbish = []
    kwargs  = dict([(kw,it) for kw,it in kwargs.items() if kw in args or rubbish.append((kw,it))])

    return fun(stream,**kwargs)
    
    
