import versions

from ftypes import FFParam, Ilist
from tprio  import *


class Block(ListWithNames):
    def __init__(self,tpr):
        self.extend([
            ("dummy", Tuple(tpr, 256, Integer) if tpr.version < 44 else None),
            ("nr",    Integer(tpr)),
            ("nra",   Integer(tpr) if tpr.version < 51 else None),
        ])
        self.extend([
            ("index", Tuple(tpr, self.nr+1, Integer)),
            ("A",     Tuple(tpr, self.nra,  Integer) if tpr.version < 51 else None )
        ]) 


class BlockA(ListWithNames):
    def __init__(self,tpr):
        self.extend([
            ("dummy", Tuple(tpr, 256, Integer) if tpr.version < 44 else None),
            ("nr",    Integer(tpr)),
            ("nra",   Integer(tpr)),
        ])
        self.extend([
            ("index", Tuple(tpr, self.nr+1, Integer)),
            ("A",     Tuple(tpr, self.nra,  Integer))
        ]) 


class Atom(ListWithNames):
    def __init__(self,tpr):

        self.position = tpr.tell()

        ngrps = versions.ngrps(tpr)

        self.extend([
            ("m",      Real(tpr)),
            ("q",      Real(tpr)),
            ("mB",     Real(tpr)),
            ("qB",     Real(tpr)),
            ("type",   Unsigned(tpr)),
            ("typeB",  Unsigned(tpr)),
            ("ptype",  Integer(tpr)),
            ("resind", Integer(tpr)),
            ("atomnumber", Integer(tpr) if tpr.version >= 52 else 0),
            ("groups",     Tuple(tpr, ngrps, Unsigned) if tpr.version < 57 else None)
            ])


class Moltype(ListWithNames):
    def __init__(self,tpr,symtab=None):
        self.position     = tpr.tell()

        ngrps = versions.ngrps(tpr)

        self.extend([
            ("nameidx",     Integer(tpr) if tpr.version >= 57 else None), # Symtab index for Moleculetype Name
            ("natoms",      Integer(tpr)),                                # Number of atoms in Moleculetype
            ("nres",        Integer(tpr)),                                # Number of residues in Moleculetype
            ("ngroupname",  Integer(tpr) if tpr.version < 57 else None),  # Number of group names in Moleculetype.
                                                                          # Kept in struct starting with version 57.
            ])

        self.extend([
            ("atoms",        Tuple(tpr, self.natoms, Atom)),
            ("atomnameidx",  Tuple(tpr, self.natoms, Integer)),
            ("atomtypeidx",  Tuple(tpr, self.natoms, Integer) if tpr.version > 20 else []),
            ("atomtypeBidx", Tuple(tpr, self.natoms, Integer) if tpr.version > 20 else []),
            ("residues",     Tuple(tpr, self.nres, (Integer,Integer,Char) if tpr.version >= 63 else (Integer,))),
            ("groupnameidx", Tuple(tpr, self.ngroupname, Integer) if tpr.version < 57 else []),
            ("groups",       Tuple(tpr, ngrps, Group) if tpr.version < 57 else []),
            ("ilists",       Ilist(tpr)  if tpr.version >= 57 else None),
            ("cgs",          Block(tpr)  if tpr.version >= 57 else None), 
            ("excls",        BlockA(tpr)),
        ])

        if symtab:
            self.setnames(symtab)

    def setnames(self,symtab):
        self.extend([
            ("name",         symtab[self.nameidx] if self.nameidx != None else None),
            ("atomnames",    [ symtab[i]    for i in self.atomnameidx  ]),
            ("atomtypes",    [ symtab[i]    for i in self.atomtypeidx  ]),
            ("atomtypesB",   [ symtab[i]    for i in self.atomtypeBidx ]),
            ("groupnames",   [ symtab[i]    for i in self.groupnameidx ]),
            ("residuenames", [ symtab[i[0]] for i in self.residues     ]),
        ])


class Molblock(ListWithNames):
    def __init__(self,tpr):
        self.position = tpr.tell()
        self.extend([
            ("type",       Integer(tpr)),
            ("nmol",       Integer(tpr)),
            ("natoms_mol", Integer(tpr)),
            ("posresA",    Group(tpr, (Real,Real,Real))),
            ("posresB",    Group(tpr, (Real,Real,Real))),
        ])


class AtomTypes(ListWithNames):
    def __init__(self,tpr):
        if tpr.version <= 25:
            return
        self.append(("atomTypeN", Integer(tpr)))
        self.extend([
            ("radius",     Tuple(tpr, self.atomTypeN, Real)),
            ("volume",     Tuple(tpr, self.atomTypeN, Real)),
            ("surftens",   Tuple(tpr, self.atomTypeN, Real)),
            ("number",     Tuple(tpr, self.atomTypeN, Integer) if tpr.version >= 40 else None),
            ("gbRadius",   Tuple(tpr, self.atomTypeN, Real)    if tpr.version >= 60 else None),
            ("S_hct",      Tuple(tpr, self.atomTypeN, Real)    if tpr.version >= 60 else None),
        ])


class Topology(ListWithNames):
    def __init__(self,tpr):
        
        ngrps = versions.ngrps(tpr)

        self.extend([
            ("symtab",       Strings(tpr)),
            ("symstridx",    Integer(tpr)),
        ])
        
        self.extend([
            ("symstr",       self.symtab[self.symstridx]),
            ("ffparam",      FFParam(tpr,self.symtab) if tpr.version >= 57 else None),
            ("moltypes",     Group(tpr, Moltype)      if tpr.version >= 57 else None),
            ("molblocks",    Group(tpr, Molblock)     if tpr.version >= 57 else None),
            ("topnatoms",    Integer(tpr)             if tpr.version >= 57 else None),
            ("atomtypes",    AtomTypes(tpr)           if tpr.version >  25 else None),

            # For earlier versions (<57), there should be one moltype, and
            # this set of ffparam/ilists should be set as attribute to that
            ("ffparam",      FFParam(tpr,self.symtab) if tpr.version <  57 else None),
            ("ilists",       Ilist(tpr)               if tpr.version <  57 else None),

            ("cmapN",        Integer(tpr)             if tpr.version >= 65 else None),
            ("cmapGridSize", Integer(tpr)             if tpr.version >= 65 else None),
        ])

        self.extend([
            ("cmap",         [ Tuple(tpr, self.cmapGridSize**2, (Real, Real, Real, Real)) for i in range(self.cmapN) ] if tpr.version >= 65 else None ),
            ("groupids",     [ Group(tpr, Integer) for i in range(ngrps) ] if tpr.version >= 57 else None ),
            ("groupnameidx", Group(tpr, Integer) if tpr.version >= 57 else None),
            ("groups",       [ Group(tpr, Unsigned) for i in range(ngrps) ] if tpr.version >= 57 else None ),
            ("cgs",          Block(tpr)  if tpr.version < 57 else None),
            ("mol",          Block(tpr)  if tpr.version < 57 else None),
            ("shake",        BlockA(tpr) if tpr.version < 51 else None),
        ])
        
        



