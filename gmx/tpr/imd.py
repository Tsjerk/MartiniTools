import versions 

from tprio import *

class InteractiveMD(versions.Test):
    def __init__(self,tpr):
        self.position = tpr.tell()
        if tpr.version < versions.InteractiveMolecularDynamics:
            self.empty = True
            return
        self.atoms = Integer(tpr) and [ Integer(tpr) for i in xrange(Integer(tpr)) ]

    def __nonzero__(self):
        return not hasattr(self,"empty")

