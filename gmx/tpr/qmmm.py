import versions
from tprio import *

class QMMM(ListWithNames):

    def __init__(self,tpr):

        self.position = tpr.tell()

        if tpr.version < 39:
            self.empty = True
            return

        self.extend([
            
            # Boolean 
            ("QMMM",         Integer(tpr)),
            
            # Scheme for QM/MM calculations
            ("QMMMscheme",    Integer(tpr)),

            # Scale factor for MM charges
            ("MMChargeScaleFactor",   Real(tpr)),

            # Number of QM groups
            ("ngQM",          Integer(tpr)),
        ])

        self.extend([

            # QM method; index to enum
            ("QMmethod",      Tuple(tpr, self.ngQM, Integer)),

            # QM basis set; index to enum
            ("QMbasis",       Tuple(tpr, self.ngQM, Integer)),

            # Charge of QM part
            ("QMcharge",      Tuple(tpr, self.ngQM, Integer)),

            # Multiplicity of QM part
            ("QMmult",        Tuple(tpr, self.ngQM, Integer)),

            # Surface hopping; boolean
            ("SH",            Tuple(tpr, self.ngQM, Integer)),

            # CAS space options            
            ("CASorbitals",   Tuple(tpr, self.ngQM, Integer)),
            ("CASelectrons",  Tuple(tpr, self.ngQM, Integer)),
            ("SAon",          Tuple(tpr, self.ngQM, Real)),
            ("SAoff",         Tuple(tpr, self.ngQM, Real)),
            ("SAsteps",       Tuple(tpr, self.ngQM, Integer)),

            # Optimization of QM subsystem
            ("bOPT",          Tuple(tpr, self.ngQM, Integer)),
            ("bTS",           Tuple(tpr, self.ngQM, Integer)),
        ])
        
    def __nonzero__(self):
        return not hasattr(self,"empty")
        
