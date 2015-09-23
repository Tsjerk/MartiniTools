import versions

class RotateGroup(versions.Test):
    def __init__(self,tpr):
        self.eType          = tpr.readInteger()
        self.bMassW         = tpr.readInteger()
        self.nat            = tpr.readInteger()
        self.ind            = [ tpr.readInteger() for i in range(self.nat) ]
        self.x_ref          = [ tpr.parse("RRR")  for i in range(self.nat) ]
        self.vec            = tpr.parse("RRR")
        self.pivot          = tpr.parse("RRR")
        self.rate           = tpr.readReal()
        self.k              = tpr.readReal()
        self.slab_dist      = tpr.readReal()
        self.min_gaussian   = tpr.readReal()
        self.eps            = tpr.readReal()
        self.eFittype       = tpr.readInteger()
        self.PotAngle_nstep = tpr.readInteger()
        self.PotAngle_step  = tpr.readReal()


class Rotate(versions.Test):
    def __init__(self,tpr):
        if tpr.version < 74 or not tpr.readInteger():
            self.empty = True
            return
        self.ngrp    = tpr.readInteger()
        self.nstrout = tpr.readInteger()
        self.nstsout = tpr.readInteger()
        self.groups  = [ RotateGroup(tpr) for i in range(self.ngrp) ]


