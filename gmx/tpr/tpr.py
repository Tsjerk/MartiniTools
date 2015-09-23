
import versions  # Version management stuff
import numpy

from inputrec import InputRecord # MD control parameters
from topology import Topology    # Molecule and Forcefield stuff
from tprio    import *           # TPR Stream reading/formatting routines

import sys 

# Check version
if sys.hexversion < 0x02070000:
    class PythonTooOldError(BaseException):
        def __init__(self,msg): self.msg = msg
        def __str__(self):      return self.msg
    raise(PythonTooOldError("\nThe TPR reader requires at least Python version 2.7, "
                            "while it's called with %s:\n%s\n"%(sys.executable,sys.version)))

class CoordSet(numpy.ndarray):
    def __new__(cls, coord, time, box):
        obj      = numpy.asarray(coord).view(cls)
        obj.box  = numpy.asarray(box)
        obj.time = time
        return obj


class TPR(ListWithNames):
    def __init__(self,filename):
        # TPRIO is an io.BytesIO derived class,
        # which has no public constructor. 
        # The TPR file is read into memory (as string)
        # and then wrapped as an io.BytesIO object to 
        # allow file-like processing. This approach 
        # avoids system calls for reading from disk.
        with open(filename) as f:
            self.tpr = f.read()

        tpr = TPRIO(self.tpr)
        tpr.init()
        version        = tpr.version

        print >>sys.stderr, "TPR file {} {}, file version {}, generation {}, precision {}".format(
            tpr.gmx, tpr.tag, tpr.version, tpr.generation, tpr.precision)
        
        if tpr.version <= 9:
            raise versions.VersionError("Cannot read TPR files older than version 10. Where did you get this? Carved on a clay tablet?")

        if tpr.version == 80:
            raise versions.VersionError("TPR version 80 is ambiguous, used for GMX 4.6 and 5.0 development versions,"
                                        " and can't be processed. Please make a new one.")

        # TPR data ===>>>

        self.extend([
            ("gmx",          String(tpr)),
            ("precision",    Integer(tpr)),
            ("version",      Integer(tpr)),
            ("tag",          String(tpr)  if 77 <= version <= 80 else None),
            ("generation",   Integer(tpr) if       version >= 26 else 0),
            ("tag",          String(tpr)  if       version >= 81 else tpr.tag),
            ("natoms",       Integer(tpr)),
            ("ngtc",         Integer(tpr) if tpr.version >= 28 else 0),
            ("dummy",        Integer(tpr) if tpr.version <  62 else 0),
            ("dummy",        Real(tpr)    if tpr.version <  62 else 0),
            ("fep_state",    Integer(tpr) if version >= 79 else 0),       # Free energy state. Will probably change in the future
            ("lambd",        Real(tpr)                           ),       # Free energy lambda
            ("hasIR",        Integer(tpr)                        ),       # Presence of input record
            ("hasTop",       Integer(tpr)                        ),       # Presence of topology
            ("hasX",         Integer(tpr)                        ),       # Presence of coordinates
            ("hasV",         Integer(tpr)                        ),       # Presence of velocities
            ("hasF",         Integer(tpr)                        ),       # Presence of forces
            ("hasBox",       Integer(tpr)                        ),       # Presence of unit cell definition
            ])

        self.extend([
            ("box",          RealMatrix(tpr) if self.hasBox else None   ),                  # Periodic boundary conditions lattice vectors (box matrix)
            ("box_rel",      RealMatrix(tpr) if version >= 51 else None ),                  # ...
            ("box_vel",      RealMatrix(tpr) if version >= 28 else None ),                  # Box vector velocities
            ("dummy",        RealMatrix(tpr) if 28 <= version <= 56 else None),             # A redundant matrix
            ("dummy",        Tuple(tpr, self.ngtc, Real) if 28 <= version < 69 else None),
            ("dummy",        Tuple(tpr, self.ngtc, Real) if 28 <= version      else None),  # Used to be Berendsen t-couple lambdas
            ("ir",           InputRecord(tpr) if self.hasIR and version < 26 else None), 
            ("top",          Topology(tpr) if self.hasTop else None),
            ("x",            Coordinates(tpr, self.natoms, self.precision) if self.hasX else None), 
            ("v",            Coordinates(tpr, self.natoms, self.precision) if self.hasV else None), 
            ("f",            Coordinates(tpr, self.natoms, self.precision) if self.hasF else None), 
            ("pbc",          Integer(tpr) if self.hasIR and tpr.version >= 53 else None),
            ("periodicMols", Integer(tpr) if self.hasIR and tpr.version >= 53 else None), 
            ("ir",           InputRecord(tpr) if self.hasIR and tpr.version >= 26 else None),
        ])

        # Dictionary with modifications: position -> modified binary string
        self.modifications = {}


    def __str__(self):        
        if not self.modifications:
            return self.tpr

        mods = self.modifications.items() + [(len(self.tpr),"")]
        mods.sort()
        positions, strings = zip(*mods)
        ends = [ start+len(s) for start,s in mods ]

        out = []
        for s,e,m in zip([0]+ends, positions, strings):
            out.append(self.tpr[s:e])
            out.append(m)

        return "".join(out)

    def get(self, param, value=None):
        par = self.search(param)
        if par is None or par is False:
            raise IndexError("Parameter not found: {}".format(param))
        return par

    def set(self, param, value=None):
        par = self.search(param)
        if par is None or par is False:
            raise IndexError("Parameter not found: {}".format(param))
        newval = par.pack(value)
        self.modifications[par.position] = newval
        return param, par, newval
