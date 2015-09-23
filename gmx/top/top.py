


import sys, random, math, re, os, itertools



# Need to extract moleculetypes, residues and atom lists


# Set the pattern for matching files to be included
includePattern = re.compile('#include "(.*)"')


# Gromacs force field directory
gmxlib = os.environ.get("GMXLIB")
if not gmxlib:
    gmxdat = os.environ.get("GMXDATA")
    if gmxdat:
        gmxlib = os.path.join(gmxdat,"gromacs","top")
    else:
        gmxlib="."



# The following function finds and follows #included files
def reciter(filename):
    # Set the directory of the filename so we know where to expect #included files
    dir = os.path.dirname(filename)

    # Iterate over the lines
    for line in open(filename):

        # Check for an #include statement; yield the line if there is none
        if line.strip().startswith("#include"):
            
            # Extract the #include filename             
            matches = re.findall(includePattern,line)

            if matches:
                fr = matches[0]

                if not os.path.exists(fr):
                    fr = os.path.join(dir,matches[0])

                if not os.path.exists(fr):
                    fr = os.path.join(gmxlib,matches[0])

                if not os.path.exists(fr):
                    yield "; " + line + " ; File not found\n"
                else:
                    for j in reciter(fr):
                        yield j
        else:
            yield line


class TOP:
    def __init__(self,other,out=None):
        
        # Process the topology file extract moleculetypes, atom lists and the molecule list
        self.molecules = []
        self.top       = []

        # Processed topology
        # This is equal to the input target topology, but with all #include statements resolved
        if out:
            out = open(out,"w")

        # List of line numbers at which to find moleculetype definitions
        mols     = []

        # Moleculetypes
        moltypes = []

        # Atoms per moleculetype
        atoms    = []

        # Gromacs topology directive
        tag      = re.compile('^ *\[ *(.*) *\]')

        # Last directive read (current)
        cur      = None

        # Set line counter
        counter = 0

        # Iterate over lines, processing #included files
        for line in reciter(other):
            if out:
                out.write(line)

            # Increment the line counter
            counter += 1

            # Add the line to the (processed) topology
            self.top.append(line)

            # Strip leading and trailing spaces
            s = line.strip()

            # Lines starting with [ indicate a directive
            if s.startswith("["):
                # Extract the directive name
                cur = re.findall(tag,s)[0].strip()
                continue

            # Conditionals :S
            # Conditionals are simply skipped             
            if s.startswith("#"):
                continue            

            # Strip comments
            s = s.split(';')[0].strip()

            # Skip empty lines
            if not s:
                continue

            if cur == "moleculetype":
                moltypes.append(s.split()[0])
                atoms.append([])
                mols.append(counter)
                continue

            if cur == "system":
                mols.append(counter)
                continue

            if cur == "atoms":
                # Comments are already skipped
                a = s.split()                
                atoms[-1].append((a[4],a[3],a[2],""))
                continue

            if cur == "molecules":
                # Each molecules entry has a moleculetype name and a number
                # The moleculetype name is added to the molecules list 
                # as many times as the number indicates. This makes it easy
                # to expand the molecules to atoms.
                m = s.split()
                for j in range(int(m[1])):
                    self.molecules.append(m[0])

        if out:
            out.close()


        # Convert moleculetypes to dictionary
        self.moleculetypes = dict(zip(moltypes,atoms))

        molecules = [(i,len(list(j))) for i,j in itertools.groupby(self.molecules)]

        # Build a full atom list
        # The chain identifier is unique for each molecule
        # The moleculetype name is added as last element
        mr = zip(self.molecules, range(len(self.molecules)))
        self.atoms = [[a,r,i,c,0,0,0,t] for t,c in mr for a,r,i,m in self.moleculetypes[t]]

        # Build a residue list
        if self.atoms:
            self.residues = [[self.atoms[0]]]        
            for i in self.atoms[1:]:
                if i[1:4] != self.residues[-1][-1][1:4]:
                    self.residues.append([])
                self.residues[-1].append(i)
        else:
            self.residues = []


