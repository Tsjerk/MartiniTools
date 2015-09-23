
from structure import Structure

import math, numpy, gzip


# Reformatting of lines in structure file                                     
pdbAtomLine = "ATOM  %5d %4s%4s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n"        
pdbBoxLine  = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n"        


# Conversion from degrees to radians
_d2r = 3.14159265358979323846264338327950288/180


def pdbOut(atom,i=1):
    insc = atom[2]>>20
    resi = atom[2]-(insc<<20)
    pdbline = "ATOM  %5i  %-3s %3s%2s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f           %1s  \n"
    return pdbline%((i,atom[0][:3],atom[1],atom[3],resi,chr(insc)) + atom[4:] + (1,40,atom[0][0])) 


def isPdbAtom(a):    
    return a.startswith("ATOM") or (options["-hetatm"] and a.startswith("HETATM")) or a.startswith("TER")


def pdbReadBox(a):
    fa, fb, fc, aa, ab, ac = [float(i) for i in a.split()[1:7]]
    fa *= 0.1 # Conversion to nm
    fb *= 0.1 # Conversion to nm
    fc *= 0.1 # Conversion to nm
    ca = math.cos(_d2r*aa)
    cb = math.cos(_d2r*ab)
    cg = math.cos(_d2r*ac)
    sg = math.sin(_d2r*ac)
    wx = fc*cb
    wy = fc*(ca-cb*cg)/sg 
    wz = math.sqrt(fc*fc - wx*wx - wy*wy)
    return numpy.array([[   fa,     0,   0], 
                        [fb*cg, fb*sg,   0], 
                        [   wx,    wy,  wz]])


def pdbBoxString(box):
    # Box vectors
    u, v, w  = box[0,:], box[1,:], box[2,:]

    # Box vector lengths
    nu,nv,nw = [math.sqrt(norm2(i)) for i in (u,v,w)]

    # Box vector angles
    alpha = nv*nw == 0 and 90 or math.acos(cos_angle(v,w))/_d2r
    beta  = nu*nw == 0 and 90 or math.acos(cos_angle(u,w))/_d2r
    gamma = nu*nv == 0 and 90 or math.acos(cos_angle(u,v))/_d2r

    # Mind the factor 10 for conversion to Angstrom
    return pdbBoxLine % (10*norm(u),10*norm(v),10*norm(w),alpha,beta,gamma)


def pdbAtom(a):
    ##01234567890123456789012345678901234567890123456789012345678901234567890123456789
    ##ATOM   2155 HH11 ARG C 203     116.140  48.800   6.280  1.00  0.00
    if a.startswith("TER"):
        return 0

    # NOTE: The 27th field of an ATOM line in the PDB definition can contain an
    #       insertion code. We shift that 20 bits and add it to the residue number
    #       to ensure that residue numbers will be unique.

    return (a[12:16],                       # Name
            a[17:20],                       # Residue name
            int(a[22:26]),#+(ord(a[26])<<20), # Residue number
            a[21],                          # Chain
            0.1*float(a[30:38]),            # X
            0.1*float(a[38:46]),            # Y
            0.1*float(a[46:54]),            # Z
            float(a[54:60]),                # Occupancy
            float(a[60:66]))                # B-factor


# Simple PDB iterator
def pdbFrameIterator(stream):  
    if type(stream) == str:
        if stream.endswith(".gz"):
            stream = gzip.open(stream)
        else:
            stream = open(stream)

    title, atoms, box = [], [], []

    for i in stream:
        if i.startswith("ENDMDL"):
            yield Structure(title="".join(title), atoms=atoms, box=box)
            title, atoms, box = [], [], None            
        elif i.startswith("TITLE"):
            title.append(i)
        elif i.startswith("CRYST1"):
            box = pdbReadBox(i)
        elif i.startswith("ATOM") or i.startswith("HETATM"):
            atoms.append(pdbAtom(i))

    if atoms:
        yield Structure(title="".join(title), atoms=atoms, box=box)


