import versions

from tprio import *

class AdResS(versions.Test):
    def __init__(self,tpr):
        self.position = tpr.tell()
        if tpr.version < 77 or not Integer(tpr):
            self.empty      = True
            return
        self.type           = Integer(tpr)
        self.const_wf       = Float(tpr, tpr.precision)
        self.ex_width       = Float(tpr, tpr.precision)
        self.hy_width       = Float(tpr, tpr.precision)
        self.icor           = Integer(tpr)
        self.site           = Integer(tpr)
        self.refs           = Vector(tpr, tpr.precision)
        self.n_tf_grps      = Integer(tpr)
        self.ex_forcecap    = Float(tpr,tpr.precision)
        self.n_energy_grps  = Integer(tpr)
        self.do_hybridpairs = Integer(tpr)
        self.tf_table_index = [ Integer(tpr) for i in range(self.n_tf_grps) ]
        self.group_explicit = [ Integer(tpr) for i in range(self.n_energy_grps) ]

    def __nonzero__(self):
        return not hasattr(self,"empty")

