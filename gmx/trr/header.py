
class TrrHeader:
    """Container for header info from Gromacs' TRR file"""

    # Header format (after tag, positions 24:84)
    # >   - big endian
    # l   - input record size (usually 0)                       ---
    # l   - energy record size (usually 0)                       |
    # l   - box size (9*4 = 36 bytes, usually always)            |
    # l   - virial size (0)                                      |
    # l   - pressure size (0)                                 frsize =
    # l   - topology size (0)                                 sum + 84
    # l   - symbol table size (0)                                |
    # l   - X array size (dim*floatSize*natoms, if present)      |
    # l   - V array size (dim*floatSize*natoms, if present)      |
    # l   - F array size (dim*floatSize*natoms, if present)     --- 
    # l   - NATOMS
    # l   - step 
    # l   - nre (number of energy terms?)
    # f/d - time (depends on TRR file precision)
    # f/d - lambda (depends on TRR file precision)

    def __init__(self, header, dim=3):
        stuff = struct.unpack('>lllllllllllll', header)
        self.inputRec    = stuff[0]
        self.energy      = stuff[1]
        self.box         = stuff[2]
        self.virial      = stuff[3]
        self.pressure    = stuff[4]
        self.topology    = stuff[5]
        self.symboltable = stuff[6]
        self.x           = stuff[7]
        self.v           = stuff[8]
        self.f           = stuff[9]
        self.atoms       = stuff[10]
        self.step        = stuff[11]
        self.nre         = stuff[12]
        self.bytesize    = sum(stuff[:10])
        self.block       = dim*self.atoms
        self.float       = self.box//(dim*dim) or (self.x or self.v or self.f)//self.block 

        # Set the dtype and the dtype with endiannes
        self.dtype = (self.float == 4 and "f") or (self.float == 8 and "d") 

        if not self.dtype:
            raise IOError("Unable to determine precision of TRR file. Float size appears to be {}".format(self.float))

        self.dtypeE = '>'+self.dtype


