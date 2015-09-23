import versions

from tprio import *

class SwapCoords(versions.Test):
    def __init__(self,tpr):
        self.position = tpr.tell()

        if tpr.version < versions.ComputationalElectrophysiology:
            self.empty = True
            return

        self.eSwapCoords = Integer(tpr)

        if not self.eSwapCoords:
            self.empty = True
            return
        
        self.natoms         = Integer(tpr)
        self.satoms         = Integer(tpr)
        
        self.natoms_split_a = Integer(tpr)
        self.massw_split_a  = Integer(tpr)

        self.natoms_split_b = Integer(tpr)
        self.massw_split_b  = Integer(tpr)

        self.nstswap        = Integer(tpr)
        self.nAverage       = Integer(tpr)
        self.threshold      = Integer(tpr)
        self.cyl0r          = Float(tpr, tpr.precision)
        self.cyl0u          = Float(tpr, tpr.precision)
        self.cyl0l          = Float(tpr, tpr.precision)
        self.cyl1r          = Float(tpr, tpr.precision)
        self.cyl1u          = Float(tpr, tpr.precision)
        self.cyl1l          = Float(tpr, tpr.precision)                 

        self.indices       = [ Integer(tpr) for i in xrange(self.natoms) ]
        self.sol_indices   = [ Integer(tpr) for i in xrange(self.satoms) ]

        self.indices_split = [[Integer(tpr) for i in xrange(j)] for j in self.natoms_split]
        
        # There are currently two compartments (GMX 5.0rc1) 
        eCompNR = 2
        self.ions = [ (Integer(tpr), Integer(tpr)) for i in xrange(eCompNR) ]
            
    def __nonzero__(self):
        return not hasattr(self, "empty")

