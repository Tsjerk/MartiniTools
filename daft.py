#!/usr/bin/env python

import sys, random, math, os, itertools, numpy

version = "TAW 20140707.2120"

pi   = math.pi
sqrt = math.sqrt

# Units are shifted over x. In most cases, this arrangement is rotated 
# with respect to the reduced (optimal) lattice. These rotations are 
# well defined and listed here (radians):
rotations = [0, pi/4, pi/6, 0, -math.atan(0.5), -math.atan(sqrt(0.12)),
             -math.atan(sqrt(0.12)), -math.atan(1./3), -math.atan(1./3)]

# The reduced lattice is well-defined for each arrangement. 
# Each lattice is given as upper triangular matrix: ux, vx, vy, wx, wy, wz
#
# NOTE:
# To avoid interactions over periodic boundaries (binding of a partner at 
# two sides) the diagonals have to be set large enough. For the case of
# two solutes this means setting basic cubic unit cell edge to
# (dA + dB + r)/sqrt(2).
#
# How this works out for higher order complexes is not really certain.
# This is accounted for by multiplying the dimension later on. The table 
# here is unaltered. 
boxes = (
    (       1,       0.5,  sqrt(3)/2,       0.5,     sqrt(3), sqrt(3./2)), #1 (0)
    ( sqrt(2),         0,    sqrt(2),  sqrt(.5),    sqrt(.5),          1), #2
    ( sqrt(3), sqrt(3)/2,        1.5, sqrt(3)/2,         0.5,          1), #3
    (       2,         0,    sqrt(3),         1,   sqrt(3)/3, sqrt(3./2)), #4
    ( sqrt(5),         0,    sqrt(5), sqrt(5)/2,   sqrt(5)/2,   sqrt(.5)), #5
    ( sqrt(7), sqrt(7)/2, sqrt(21)/4, sqrt(7)/2, sqrt(21)/12, sqrt(3./2)), #6 (7)
    ( sqrt(7), sqrt(7)/2, sqrt(21)/4, sqrt(7)/2, sqrt(21)/12, sqrt(3./2)), #7
    (sqrt(10),3*math.cos(pi+rotations[-1]),3*math.sin(pi+rotations[-1]),
     sqrt(4.5)*math.cos(pi/2+rotations[-1]), sqrt(4.5)*math.sin(pi/2+rotations[-1]),sqrt(.5)),
    (sqrt(10),3*math.cos(pi+rotations[-1]),3*math.sin(pi+rotations[-1]),
     sqrt(4.5)*math.cos(pi/2+rotations[-1]), sqrt(4.5)*math.sin(pi/2+rotations[-1]),sqrt(.5))
    )

c = math.cos(pi/6)

positions = [
    [],                                                                                                            # 0
    [ (0.0, 0.0) ],                                                                                                # 1
    [ (0.0, 0.0), (1.0, 0.0) ],                                                                                    # 2
    [ (0.0, 0.0), (1.0, 0.0), (0.5,   c) ],                                                                        # 3
    [ (0.0, 0.0), (1.0, 0.0), (0.5,   c), (1.5,   c) ],                                                            # 4
    [ (1.0, 0.0), (2.0, 0.0), (1.0, 1.0), (1.0,-1.0), (0.0, 0.0) ],                                                # 5
    [ (0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (0.5,   c), (1.5,   c), (0.5,  -c), (1.5,  -c)],                         # 6
    [ (0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (0.5,   c), (1.5,   c), (0.5,  -c), (1.5,  -c)],                         # 7
    [ (1.0, 1.0), (2.0, 1.0), (0.0, 2.0), (1.0, 2.0), (2.0, 2.0), (0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (0.0, 1.0)], # 8
    [ (1.0, 1.0), (2.0, 1.0), (0.0, 2.0), (1.0, 2.0), (2.0, 2.0), (0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (0.0, 1.0)], # 9
    ]

vectors = [
    ((1.0, 0.0, 0,0), ( 0.5,   c, 0.0), (0.5, c/3, c*math.sqrt(8./9))), # 0 Rhombic dodecahedron with hexagonal xy-plane
    ((1.0, 0.0, 0,0), ( 0.5,   c, 0.0), (0.5, c/3, c*math.sqrt(8./9))), # 1 Rhombic dodecahedron with hexagonal xy-plane
    ((1.0,-1.0, 0.0), ( 1.0, 1.0, 0.0), (1.0, 0.0, 1.0)),               # 2
    ((1.5,  -c, 0.0), ( 0.0, 2*c, 0.0), (1.0, 0.0, 1.0)),               # 3
    ((2.0, 0.0, 0.0), ( 0.0, 2*c, 0.0), (1.0,   c, 1.0)),               # 4
    ((2.0,-1.0, 0.0), ( 1.0, 2.0, 0.0), (1.5, 0.5, 1.0)),               # 5
    ((2.5,  -c, 0.0), (-0.5, 3*c, 0.0), (1.5, c/3, 1.0)),               # 6
    ((2.5,  -c, 0.0), (-0.5, 3*c, 0.0), (1.5, c/3, 1.0)),               # 7
    ((3.0,-1.0, 0.0), ( 0.0, 3.0, 0.0), (2.0,1./3, 1.0)),               # 8
    ((3.0,-1.0, 0.0), ( 0.0, 3.0, 0.0), (2.0,1./3, 1.0)),               # 8
    ]

vectors = [
    ((1.0, 0.0, 0,0), ( 0.5,   c, 0.0), (0.5, c/3, c*math.sqrt(8./9))), # 0 Rhombic dodecahedron with hexagonal xy-plane
    ((1.0, 0.0, 0,0), ( 0.5,   c, 0.0), (0.5, c/3, c*math.sqrt(8./9))), # 1 Rhombic dodecahedron with hexagonal xy-plane
    ((1.0,-1.0, 0.0), ( 1.0, 1.0, 0.0), (1.0, 0.0, 1.0)),               # 2
    ((1.5,  -c, 0.0), ( 0.0, 2*c, 0.0), (0.0, 0.0, 1.0)),               # 3
    ((2.0, 0.0, 0.0), ( 0.0, 2*c, 0.0), (0.0, 0.0, 1.0)),               # 4
    ((2.0,-1.0, 0.0), ( 1.0, 2.0, 0.0), (0.0, 0.0, 1.0)),               # 5
    ((2.5,  -c, 0.0), (-0.5, 3*c, 0.0), (0.0, 0.0, 1.0)),               # 6
    ((2.5,  -c, 0.0), (-0.5, 3*c, 0.0), (0.0, 0.0, 1.0)),               # 7
    ((3.0,-1.0, 0.0), ( 0.0, 3.0, 0.0), (0.0, 0.0, 1.0)),               # 8
    ((3.0,-1.0, 0.0), ( 0.0, 3.0, 0.0), (0.0, 0.0, 1.0)),               # 8
    ]

###

class Scheme:
    def __init__(self,pbc,grid):
        self.pbc  = numpy.array(pbc).reshape((3,3))
        self.grid = numpy.concatenate((grid,numpy.zeros((1,grid.shape[1]))))
        
    def size(self,r,z=None):
        pbc, grid = r*self.pbc, r*self.grid
        if z:
            pbc[2,2] = z
        self.grid[2,:] = pbc[2,2]/2
        return pbc, grid


            

def rot(v):
    '''Determine rotation matrix from vector, aligning vector with x'''
    n = math.sqrt(v[0]**2+v[1]**2)
    c,s = v[0]/n, v[1]/n
    return ((c,s,0),(-s,c,0),(0,0,1))


def vvadd(a,b):    
    if type(b) in (int,float):
        return [i+b for i in a]
    return [i+j for i,j in zip(a,b)]


def norm2(a):
    return sum([i*i for i in a])


def norm(a):
    return math.sqrt(norm2(a))


def iprod(a,b):
    return sum([i*j for i,j in zip(a,b)])


def mvmul(A,b):
    return [iprod(a,b) for a in A]


def isPDBAtom(l):
    return l.startswith("ATOM") or l.startswith("HETATM")


def pdbAtom(a):
    ##01234567890123456789012345678901234567890123456789012345678901234567890123456789
    ##ATOM   2155 HH11 ARG C 203     116.140  48.800   6.280  1.00  0.00
    ## ===>  cc atom name,cc   res name, c    res id, chain,
    return (str(a[12:16]),str(a[17:20]),int(a[22:26]),a[21],
            #           x,                 y,                 z       
            float(a[30:38])/10,float(a[38:46])/10,float(a[46:54])/10)


d2r = 3.14159265358979323846264338327950288/180
pdbBoxLine  = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1"        


def cos_angle(a,b):
    p = sum([i*j for i,j in zip(a,b)])
    q = math.sqrt(sum([i*i for i in a])*sum([j*j for j in b]))
    return min(max(-1,p/q),1)


def pdbBoxString(b):
    # Box vectors
    u, v, w  = (b[0],b[3],b[4]), (b[5],b[1],b[6]), (b[7],b[8],b[2])

    # Box vector lengths
    nu,nv,nw = [math.sqrt(norm2(i)) for i in (u,v,w)]

    # Box vector angles
    alpha = nv*nw == 0 and 90 or math.acos(cos_angle(v,w))/d2r
    beta  = nu*nw == 0 and 90 or math.acos(cos_angle(u,w))/d2r
    gamma = nu*nv == 0 and 90 or math.acos(cos_angle(u,v))/d2r

    return pdbBoxLine % (10*norm(u),10*norm(v),10*norm(w),alpha,beta,gamma)


def groAtom(a):
    #012345678901234567890123456789012345678901234567890
    #    1PRN      N    1   4.168  11.132   5.291
    ## ===>    atom name,     res name,     res id, chain,            x,              y,              z       
    return (str(a[10:15]), str(a[5:10]), int(a[:5]), " ", float(a[20:28]),float(a[28:36]),float(a[36:44]))


class Structure:
    def __init__(self,other):
        if type(other) == str:
            t = other.find("=") 
            self.name = None
            if t > -1:
                self.name = other[:t]
                other = other[t+1:]            
            f = open(other)
            lines = f.readlines()
            f.close()
        else:
            lines = other
            self.name = None
        # Try extracting PDB atom/hetatm definitions
        rest   = []
        self.atoms  = [pdbAtom(i) for i in lines if isPDBAtom(i) or rest.append(i)]
        if not self.atoms:             
            # This should be a GRO file
            self.atoms = [groAtom(i) for i in lines[2:-1]]
        self.x, self.y, self.z = zip(*self.atoms)[4:7]
        # Center XY
        mx,my,mz    = [sum(i)/len(i) for i in (self.x,self.y,self.z)]
        self.x0     = [ i-mx for i in self.x ]
        self.y0     = [ i-my for i in self.y ]
        self.z0     = [ i-mz for i in self.z ]
        # Calculate diameter
        self.d      = 2*math.sqrt(max([i*i+j*j for i,j in zip(self.x0,self.y0)]))
        self.d3     = 2*math.sqrt(max([i*i+j*j+k*k for i,j,k in zip(self.x0,self.y0,self.z0)]))
        self.zrange = (min(self.z0),max(self.z0)) 

    def __len__(self):
        return len(self.x)

    def randrotate3D(self):
        # Random rotation in three dimensions, using quaternions
        a,  b,  c       = random.random(), 2*math.pi*random.random(), 2*math.pi*random.random()
        s,  t           = math.sqrt(1-a), math.sqrt(a)
        # The scalar part of the quaternion is multiplied by two to save an operation in the for loop
        qw, qx, qy, qz  = 2*s*math.sin(b), s*math.cos(b), t*math.sin(c), t*math.cos(c)
        qq              = 0.25*qw*qw-qx*qx-qy*qy-qz*qz         

        out = []
        for x,y,z,a in zip(self.x0,self.y0,self.z0,self.atoms):
            qp = 2*(qx*x + qy*y + qz*z)
            out.append((
                u*D + qp*qx + qq*x + qw*(qy*z-qz*y),
                v*D + qp*qy + qq*y + qw*(qz*x-qx*z),
                qp*qz + qq*z + qw*(qx*y-qy*x), a
            ))
        return out
        

# Very simple option class
class Option:
    def __init__(self,func=str,num=1,default=None,description=""):
        self.func        = func
        self.num         = num
        self.value       = default
        self.description = description
    def __nonzero__(self): 
        if self.func == bool:
            return self.value != None
        return bool(self.value)
    def __str__(self):
        return self.value and str(self.value) or ""
    def setvalue(self,v):
        if len(v) == 1:
            self.value = self.func(v[0])
        else:
            self.value = [ self.func(i) for i in v ]


# 
items = []
combi = []

# Description
desc = ""


# Option list
options = [
#   option           type         number   default description
    ("-f",    Option(items.append, 1,        None, "Input GRO or PDB file 1: Protein")),
    ("-o",    Option(str,          1,        None, "Output PDB file")),
    ("-ndx",  Option(str,          1,        None, "Index file with energygroups")),
    ("-d",    Option(float,        1,           1, "Distance between units")),
    ("-z",    Option(float,        1,           1, "z-distance")),
    ("-Z",    Option(float,        1,           1, "z-distance")),
    ("-g",    Option(float,        1,           0, "Grid spacing for placement")),
    ("-n",    Option(int,          1,           1, "Number of structures")),
    ("-name", Option(str,          1,        None, "Name of assay")),
    ("-c",    Option(combi.append, 1,        None, "Combinations")),
    ("-D",    Option(bool,         0,        None, "Set distance over diagonal, avoiding direct PBC interactions")),
    ("-F",    Option(str,          1,        None, "Parameter file")),
    ("-3d",   Option(bool,         0,        None, "Rotations in three dimensions")),
    ("-ads",  Option(float,        1,        None, "Membrane adsorption at distance specified")),
    ("-mem",  Option(float,        1,         5.0, "Membrane thickness for adsorption")),
    ]

args = sys.argv[1:]

if '-h' in args or '--help' in args:
    print "\n",__file__
    print desc or "\nSomeone ought to write a description for this script...\n"
    for thing in options:
        print type(thing) != str and "%10s  %s"%(thing[0],thing[1].description) or thing
    print
    sys.exit()


# Convert the option list to a dictionary, discarding all comments
options = dict([i for i in options if not type(i) == str])


# Process the command line
while args:
    ar = args.pop(0)
    options[ar].setvalue([args.pop(0) for i in range(options[ar].num)])


# If there is a parameter file, process that too
if options["-F"]:
    f = open(options["-F"].value)
    stuff = f.read().split()
    f.close()
    while stuff:
        ar = stuff.pop(0)
        options[ar].setvalue([stuff.pop(0) for i in range(options[ar].num)])
        

if not options["-o"]:
    options["-o"].setvalue([options["-name"].value+".pdb"])
if not options["-ndx"]:
    options["-ndx"].setvalue([options["-name"].value+".ndx"])


# Convert all structures
items = [Structure(i) for i in items]


# Find maximum XY diameter and determine minimal cell size 
d = options["-3d"] and [i.d3 for i in items] or [i.d for i in items]
if len(d) > 1:
    d.sort()
    dmin = (d.pop()+d.pop())/2
elif d:
    dmin = d.pop()
else:
    dmin = 0

# Get the total z-range to set the distance over PBC
zrange = zip(*[i.zrange for i in items])
if not options["-ads"].value is None:
    zdim = options["-mem"].value + max(zrange[1]) - min(zrange[0]) + 2*options["-ads"].value
else:
    zdim = max(zrange[1]) - min(zrange[0]) + options["-z"].value
halfz  = zdim / 2

# Set the distance between grid positions based on grid-spacing or dmin and distance
if options["-g"]:
    if options["-g"].value < dmin:
        print "Grid spacing (%f) smaller than sum of largest radii. This will cause overlaps."
        raise StupidOxError
    D = options["-g"].value
else:
    if options["-D"]:
        D = math.sqrt(2)*(dmin + options["-d"].value)
    else:
        D = dmin + options["-d"].value


# Set grid positions in XY plane
cospi_6  = 0.8660254037844387 # math.cos(math.pi/6)
tanpi_6  = 0.5773502691896257 # math.tan(math.pi/6)
 

# Make combinations

if not combi:
    if len(items) == 1:
        combi = ["1"]
    else:
        combi = ["%d!"%len(items)]

rli   = range(len(items))
stuff = set()
for i in combi:
    if i.isdigit():
        # Make only equal type combinations
        n = int(i)
        stuff.update([tuple([i for j in range(n)]) for i in rli])
    elif i[-1] == "!":
        # Make all combinations of int(n[:-1]) items        
        # This requires excluding equal type combinations
        n = int(i[:-1])
        equalItems   = [tuple([i for j in range(n)]) for i in rli]
        unequalItems = set(itertools.combinations(rli,n))
        stuff.update(unequalItems.difference(equalItems))        
    elif i.startswith("@"):
        # Item indexing, comma separated list
        stuff.add(tuple([ int(j.replace("@",""))-1 for j in i.split(",") ]))
        

# Set the positions and vectors

_rotations = [ rot(i[0]) for i in vectors ]
_positions = [ [(r[0][0]*i[0]+r[0][1]*i[1],r[1][0]*i[0]+r[1][1]*i[1]) for i in p] for r,p in zip(_rotations, positions) ]
_vectors   = [ [mvmul(r,u) for u in v] for r,v in zip(_rotations, vectors) ]


# Write a list of structures with IDs and tags


combi = list(stuff)
combi.sort()
for combination in combi:
    tags = [items[j].name or str(j) for j in combination]
    assayDir = "-".join(tags)

    # Make a directory for this combination
    if not os.path.exists(assayDir):
        os.mkdir(assayDir)
    if options["-name"]:
        assayDir = os.path.join(assayDir,options["-name"].value)
        if not os.path.exists(assayDir):
            os.mkdir(assayDir)

    # Get the (unscaled) positions for the number of items
    n = len(combination)
    if n > 9:
        print "The maximum number of chains for DAFT analysis is 9."
        raise ValueError
    if not n:
        n = 0
 
    positions = _positions[n]

    # Get the box for the number of items
    ux,uy,uz,vx,vy,vz,wx,wy,wz = [j*D for i in _vectors[n] for j in i]
    if options["-3d"] and n == 2:
        # In the case of two-component 3D daft, the Z-dimension 
        # is equal to the distance in the XY plane
        box   = (ux,vy,ux/math.sqrt(2),0,0,vx,0,wx,wy)
    else:
        box   = (ux,vy,zdim,0,0,vx,0,wx,wy)


    ## Output

    ndx = open(os.path.join(assayDir,options["-ndx"].value),"w")
    c = 65
    i = 1
    for k,t in zip(combination,tags):
        ndx.write("[ %s_%s ]\n"%(chr(c),t))
        ndx.write("\n".join([str(i+j) for j in range(len(items[k]))])+"\n")
        i += j+1
        c += 1
    ndx.close()

    if options["-o"].value:
        base = options["-o"].value[:options["-o"].value.rfind(".")]
    else:
        base = options["-name"].value
    base = os.path.join(assayDir,base)+"-%04d"
    fnm  = options["-o"] and options["-o"].value or base+".pdb"

    # Skip existing directories; allow extending existing runs
    num     = 1
    while os.path.exists(base%num):
        num += 1

    # Print the fnm, the item numbers and the starting number of the directory
    print "@ %s %d %d :"%(os.path.join(base,fnm),num,num+options["-n"].value), " ".join([str(j+1) for j in combination])   

    groline = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"
    pdbline = "ATOM  %5i %-3s %3s%2s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n"

    for i in xrange(options["-n"].value):
        id  = 1
        p   = random.sample(positions,n)

        itemDir = base%num
        os.mkdir(itemDir)
        
        f = open(os.path.join(itemDir,fnm),"w")
        f.write("MODEL %8d\n"%num)
        f.write(pdbBoxString(box)+"\n")

        for j in [items[k] for k in combination]:
            u,v = p.pop()

            # Random rotation
            if options["-3d"]:
                
                for x,y,z,a in j.randrotate3D():
                    if not options["-ads"].value is None:
                        z += halfz
                    f.write(pdbline%(id,a[0],a[1],a[3],a[2]," ",10*x,10*y,10*z,1,0))
                    id += 1

            elif not options["-ads"].value is None:
                # Random rotation around X (parallel to membrane)
                r   = random.random()*math.pi*2
                c,s = math.cos(r), math.sin(r)        
       
                for x,y,z,a in zip(j.x0,j.y0,j.z,j.atoms):
                    #     position    rotation
                    yn =    u*D   +   c*y-s*z
                    zn =    v*D   +   s*y+c*z

                    # Print the line; mind the conversion from nanometer to angstrom
                    f.write(pdbline%(id,a[0],a[1],a[3],a[2]," ",10*x,10*yn,10*(zn+halfz),1,0))
                    id += 1

            else:
                # Random rotation in XY-plane
                r   = random.random()*math.pi*2
                c,s = math.cos(r), math.sin(r)        
       
                for x,y,z,a in zip(j.x0,j.y0,j.z,j.atoms):
                    #     position    rotation
                    xn =    u*D   +   c*x-s*y
                    yn =    v*D   +   s*x+c*y

                    # Print the line; mind the conversion from nanometer to angstrom
                    f.write(pdbline%(id,a[0],a[1],a[3],a[2]," ",10*xn,10*yn,10*z,1,0))
                    id += 1

            f.write("TER\nENDMDL\n")
        f.close()
        num += 1









