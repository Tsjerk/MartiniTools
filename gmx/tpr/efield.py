import versions

from tprio import *

class Efield1(versions.Test):
    def __init__(self,tpr):
        self.position = tpr.tell()
        nx, nt        = Integer(tpr), Integer(tpr)
        if not nx:
            self.empty = True
        self.ex_a   = Tuple(tpr, nx, Real) 
        self.ex_phi = Tuple(tpr, nx, Real) 
        self.et_a   = Tuple(tpr, nx, Real) 
        self.et_phi = Tuple(tpr, nx, Real) 

    def __nonzero__(self):
        return not hasattr(self,"empty")


class Efield(versions.Test):
    def __init__(self,tpr):
        self.position = tpr.tell()
        self.x = Efield1(tpr)
        self.y = Efield1(tpr)
        self.z = Efield1(tpr)

    def __nonzero__(self):
        return not hasattr(self,"empty")

