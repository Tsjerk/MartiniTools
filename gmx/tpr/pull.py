from tprio import *

class PullGroup(ListWithNames):
    def __init__(self,tpr):
        self.position = tpr.tell()
        self.extend([
            ("atoms",   Group(tpr, Integer)),
            ("weights", Group(tpr, Real)),
            ("pbcatom", Integer(tpr)),
            ("vec",     RealVector(tpr)),
            ("init",    RealVector(tpr)),
            ("rate",    Real(tpr)),
            ("k",       Real(tpr)),
            ("kB",      Real(tpr) if tpr.version >= 56 else None),
         ])

class Pull(ListWithNames):
    def __init__(self,tpr):
        self.position = tpr.tell()

        self.append(("ePull", Integer(tpr) if tpr.version >= 48 else 0))

        if not self.ePull:
            self.empty = True
            return
            
        self.append(("ngrp", Integer(tpr)))

        self.extend([
            ("eGeom",      Integer(tpr)),
            ("dim",        Tuple(tpr, 3, Integer)),
            ("cyl_r1",     Real(tpr)),
            ("cyl_r0",     Real(tpr)),
            ("constr_tol", Real(tpr)),
            ("nstxout",    Integer(tpr)),
            ("nstfout",    Integer(tpr)),
            ("groups",     Group(tpr, self.ngrp+1, PullGroup)),
            ])
