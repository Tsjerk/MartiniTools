
# Latest version (just make sure that it is at least as high as the latest)
last = 1024

# Removed obsolete parameters
RemoveObsoleteParameters1 = 100

# interactive molecular dynamics (IMD) 
InteractiveMolecularDynamics = 99                        

# potentials for supporting coarse-grained force fields 
RestrictedBendingAndCombinedAngleTorsionPotentials = 98 

# change ld_seed from int to gmx_int64_t 
Use64BitRandomSeed = 97                                 

# support for ion/water position swaps (computational electrophysiology) 
ComputationalElectrophysiology = 96                     


# Groups atom can belong to
# The upper line lists the names used in the C code. It makes more sense to have the names as used in the MDP file.
tprGroups = ("TC",      "ENER",       "ACC",      "FREEZE",     "User1",      "User2",      "VCM",       "XTC",      "ORFIT",        "QMMM")
tprGroups = ("tc_grps", "energygrps", "acc_grps", "freezegrps", "user1_grps", "user2_grps", "comm_grps", "xtc_grps", "orire_fitgrp", "QMMM_grps")    


def ngrps(tpr):
    if tpr.version < 23:
        return 8
    elif tpr.version < 39:
        return 9
    else:
        return len(tprGroups)


class VersionError(Exception):
    def __init__(self,msg):
        self.msg = msg
    def __str__(self):
        return self.msg


class Test:
    def test(self,t,version=None):
        version = version or self.version
        if type(t) == int:
            return version < -t if t < 0 else version >= t
        elif type(t) == str:
            return getattr(self,t)
        else:
            return all([self.test(i,version) for i in t])      


