import versions

class Wall(versions.Test):
    def __init__(self,tpr):
        self.nwall      = tpr.readInteger()
        if not self.nwall:
            self.empty = True
        self.type       = tpr.readInteger()
        self.r_linpot   = tpr.readReal() if tpr.version >= 50 else -1
        self.atomtype   = (tpr.readInteger(), tpr.readInteger())
        self.density    = (tpr.readReal(), tpr.readReal())
        self.ewald_zfac = tpr.readReal()

