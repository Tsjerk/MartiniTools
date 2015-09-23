
import versions, struct

from tprio import *

# Function types, parameter specification and version of introduction
# The ones marked with '##' are types from the topology.
ftypes = [
    ('F_BONDS',             'RRRR',      0), ## 
    ('F_G96BONDS',          'RRRR',      0), ## 
    ('F_MORSE',             'RRRRRR',    0), ##  
    ('F_CUBICBONDS',        'RRR',      20), ## 
    ('F_CONNBONDS',         '',         20), ## 
    ('F_HARMONIC',          'RRRR',     20), ## 
    ('F_FENEBONDS',         'RR',       34), ## 
    ('F_TABBONDS',          'RIR',      43), ## 
    ('F_TABBONDSNC',        'RIR',      43), ## 
    ('F_RESTRBONDS',        'RRRRRRRR', 70), ## 
    ('F_ANGLES',            'RRRR',      0), ## 
    ('F_G96ANGLES',         'RRRR',      0), ## 
    ('F_RESTRANGLES',       'RR',       versions.RestrictedBendingAndCombinedAngleTorsionPotentials), ##
    ('F_LINEAR_ANGLES',     'RRRR',     76), ## 
    ('F_CROSS_BOND_BONDS',  'RRR',      30), ## 
    ('F_CROSS_BOND_ANGLES', 'RRRR',     30), ## 
    ('F_UREY_BRADLEY',      'RRRRRRRR', 30), ## 
    ('F_QUARTIC_ANGLES',    'RRRRRR',   34), ## 
    ('F_TABANGLES',         'RIR',      43), ## 
    ('F_PDIHS',             'RRRRI',     0), ## 
    ('F_RBDIHS',            12*'R',      0), ## 
    ('F_RESTRDIHS',         'RR',       versions.RestrictedBendingAndCombinedAngleTorsionPotentials), ##
    ('F_CBTDIHS',           'RRRRRR',   versions.RestrictedBendingAndCombinedAngleTorsionPotentials), ##
    ('F_FOURDIHS',          12*'R',     26), ## 
    ('F_IDIHS',             'RRRR',      0), ## 
    ('F_PIDIHS',            'RRRRI',    26), ## 
    ('F_TABDIHS',           'RIR',      43), ## 
    ('F_CMAP',              'II',       65), ## 
    ('F_GB12',              'RRRRR',    60), ## 
    ('F_GB13',              'RRRRR',    61), ## 
    ('F_GB14',              'RRRRR',    61), ## 
    ('F_GBPOL',             '',         72), 
    ('F_NPSOLVATION',       '',         72), 
    ('F_LJ14',              'RRRR',      0), ## 
    ('F_COUL14',            'RRRR',      0), 
    ('F_LJC14_Q',           'RRRRR',    41), ## 
    ('F_LJC_PAIRS_NB',      'RRRR',     41), ## 
    ('F_LJ',                'RR',        0), ## 
    ('F_BHAM',              'RRR',       0),
    ('F_LJ_LR',             '',          0), 
    ('F_BHAM_LR',           '',         32),
    ('F_DISPCORR',          '',          0),
    ('F_COUL_SR',           '',          0),
    ('F_COUL_LR',           '',          0),
    ('F_RF_EXCL',           '',         32),
    ('F_COUL_RECIP',        '',         32),
    ('F_LJ_RECIP',          '',         93),
    ('F_DPD',               '',         46),
    ('F_POLARIZATION',      'R',        30), ##
    ('F_WATER_POL',         'RRRRRR',    0), ## NOTE: Versions prior to 31 with water polarization are not supported
    ('F_THOLE_POL',         'RRRR',     36), ##
    ('F_ANHARM_POL',        'RRR',      76), ##
    ('F_POSRES',            'VVVV',      0), ##
    ('F_FBPOSRES',          'IVRR',     90), ##
    ('F_DISRES',            'IIRRRR',    0), ##
    ('F_DISRESVIOL',        '',         22),
    ('F_ORIRES',            'IIIRRR',   22), ##
    ('F_ORIRESDEV',         '',         22),
    ('F_ANGRES',            'RRRRI',     0), ##
    ('F_ANGRESZ',           'RRRRI',     0), ##
    ('F_DIHRES',            'IIRRR',    26), ##
    ('F_DIHRESVIOL',        '',         26),
    ('F_CONSTR',            'RR',        0), ##
    ('F_CONSTRNC',          'RR',        0), ##
    ('F_SETTLE',            'RR',        0), ##
    ('F_VSITE2',            'R',         0), ##
    ('F_VSITE3',            'RR',        0), ## 
    ('F_VSITE3FD',          'RR',        0), ##
    ('F_VSITE3FAD',         'RR',        0), ##
    ('F_VSITE3OUT',         'RRR',       0), ##
    ('F_VSITE4FD',          'RRR',       0), ##
    ('F_VSITE4FDN',         'RRR',      49), ##
    ('F_VSITEN',            'IR',       50), ##
    ('F_COM_PULL',          '',         46),
    ('F_EQM',               '',         20),
    ('F_EPOT',              '',          0),
    ('F_EKIN',              '',          0),
    ('F_ETOT',              '',          0),
    ('F_ECONSERVED',        '',         46),
    ('F_TEMP',              '',          0),
    ('F_VTEMP_NOLONGERUSED','',         69),
    ('F_PDISPCORR',         '',         66),
    ('F_PRES',              '',          0),
    ('F_DVDL_CONSTR',       '',         54),
    ('F_DVDL',              '',          0),
    ('F_DKDL',              '',          0),
    ('F_DVDL_COUL',         '',         79),
    ('F_DVDL_VDW',          '',         79),
    ('F_DVDL_BONDED',       '',         79),
    ('F_DVDL_RESTRAINT',    '',         79),
    ('F_DVDL_TEMPERATURE',  '',         79),
    ]


def versionFtypes(version=None):
    '''Return ftypes as they were available in the version specified''' 

    types = [ (i,f,p,v) for i,(f,p,v) in enumerate(ftypes) ]

    if version:
        for i,f,p,v in types:
            # No B state for Morse potentials
            if f == 'F_MORSE' and version < 79:
                types[i] = (i,f,'RRR',v)

            # No B state for Urey-Bradley potentials
            if f == 'F_UREY_BRADLEY' and version < 79:
                types[i] = (i,f,'RRRR',v)

            # Something stored wrong
            if f in ('F_ANGRES','F_ANGRESZ') and version < 42:
                types[i] = (i,f,'RRRR',v)

            # No B state for position restraints
            if f == 'F_POSRES' and version < 27:
                types[i] = (i,f,'VV',v)

            # No B state for RB dihedrals
            if f == 'F_RBDIHS' and version < 25:
                types[i] = (i,f,'RRRRRR',v)

            if f in ('F_GB12','F_GB13','F_GB14') and version < 68:
                types[i] = (i,f,'RRRRRRRRR',v),            

            if f == 'F_DIHRES' and version >=82:
                types[i] = (i,f,'RRRRRR',v)


    # Remove the functiontypes that were introduced later.
    # For bookkeeping the ID according to the newer version (here 5.0rc1) is kept.
    return [ i for i in types if version == None or i[3] <= version ]


class FFParam(ListWithNames):
    def __init__(self,tpr,symtab=None):
        self.position  = tpr.tell()
        fTypes         = versionFtypes(tpr.version)        
        self.atomnr    = Integer(tpr)
        self.dummy     = Integer(tpr) if tpr.version < 57 else 0
        self.ntypes    = Integer(tpr) 
        self.functypes = [ fTypes[i] for i in struct.unpack(">"+self.ntypes*"l",tpr.read(4*self.ntypes)) ]        
        self.reppow    = Float(tpr, "double") if tpr.version >= 66 else 12.0
        self.fudgeQQ   = Real(tpr) if tpr.version >= 57 else 0
        self.params    = [ tpr.parse(i[2]) for i in self.functypes ] # Index 2 is function parameter specification 
        self.fudgeQQ   = Real(tpr) if 54 <= tpr.version < 57 else self.fudgeQQ


class Ilist(ListWithNames):
    def __init__(self,tpr):
        self.position = tpr.tell()
        for i in versionFtypes(tpr.version):
            self.append(("dummy", Tuple(tpr, 256, Integer) if tpr.version < 44 else None))
            self.append(("ilist", Group(tpr, Integer)))
 
