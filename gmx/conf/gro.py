
from structure import Structure

import numpy, gzip

groline = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"                                    


def groBoxRead(a):    
    b = [float(i) for i in a.split()] + 6*[0]  # Padding for rectangular boxes

    return numpy.array([[b[0],b[3],b[4]],
                        [b[5],b[1],b[6]],
                        [b[7],b[8],b[2]]])


def groAtom(a,prec=None):
    #012345678901234567890123456789012345678901234567890
    #    1PRN      N    1   4.168  11.132   5.291

    # In PDB files, there might by an insertion code. To handle this, we internally add
    # constant to all resids. To be consistent, we have to do the same for gro files.
    # 32 equal ord(' '), eg an empty insertion code
    constant = 32<<20

    if not prec:
        prec = len(a.split(".")[-2])+1

    x,y,z = 20+prec, 20+2*prec, 20+3*prec

    return (a[10:15],            # Atom name
            a[5:10],             # Residue name
            int(a[:5]),#+constant, # Residue number 
            " ",                 # No chain information in GRO files
            float(a[20:x]),      # X
            float(a[x:y]),       # Y
            float(a[y:z]))       # Z



# Simple GRO iterator
def groFrameIterator(stream):
    """Read a GRO file stream frame by frame"""

    if type(stream) == str:
        if stream.endswith(".gz"):
            stream = gzip.open(stream)
        else:
            stream = open(stream)

    # For efficiency, we need to know the precision
    precision = None

    while True:
        try:
            title = stream.next()
        except StopIteration:
            break

        natoms = stream.next().strip()
        if not natoms:
            break

        natoms = int(natoms)
        
        if not precision:
            a = stream.next()
            prec = len(a.split(".")[-2])+1
            atoms = [groAtom(a,precision)]
            atoms.extend([groAtom(stream.next(),precision)  for i in range(natoms-1)]) 
        else:
            atoms  = [groAtom(stream.next(),precision)  for i in range(natoms)] 
        box    = groBoxRead(stream.next())

        yield Structure(title=title, atoms=atoms, box=box)



