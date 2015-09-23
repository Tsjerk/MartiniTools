#!/usr/bin/env python

import sys,math

version="dev-v01"
authors=["Helgi I. Ingolfsson", "Tsjerk A. Wassenaar"]


# Very simple option class
class Option:
    def __init__(self,func=str,num=1,default=None,description=""):
        self.func        = func
        self.num         = num
        self.value       = default
        self.description = description
    def __nonzero__(self): 
        if self.func == bool:
            return self.value != False
        return bool(self.value)
    def __str__(self):
        return self.value and str(self.value) or ""
    def setvalue(self,v):
        if len(v) == 1:
            self.value = self.func(v[0])
        else:
            self.value = [ self.func(i) for i in v ]



# Description
desc = """
This scripts creates a customized Martini lipid topology based on the head, linker and
tail specification strings provided. The topologyfollows the standard Martini 2.0 lipid
definitions as described in:
 -S.J. Marrink, A.H. de Vries, A.E. Mark.
  Coarse grained model for semi-quantitative lipid simulations.
  JPC-B, 108:750-760, 2004.
- S.J. Marrink, H.J. Risselada, S. Yefimov, D.P. Tieleman, A.H. de Vries.
  The MARTINI force field: coarse grained model for biomolecular simulations.
  J Phys Chem B, 111:7812-7824, 2007.
 -S.J. Marrink, A.H. de Vries, T.A. Harroun, J. Katsaras, S.R. Wassall.
  Cholesterol shows preference for the interior of polyunsaturated lipid membranes.
  JACS, 130:10-11, 2008.
 -H.J. Risselada, S.J. Marrink.
  The molecular face of lipid rafts in model membranes.
  PNAS, 105:17367-17372, 2008.
 -S. Baoukina, L. Monticelli, H.J. Risselada, S.J. Marrink, D.P. Tieleman.
  The molecular mechanism of lipid monolayer collapse.
  PNAS, 105:10803-10808, 2008.
??? Add insane lipid paper when readdy

WARNING:
  This script can generate topologys for numerous lipids many of wich are unrealistic
  and untested, please use with discression. 


The lipid descriptions supported are as follows:

Heads (-he): 
  Please provide a list of lipid head beads. The left most bead will be on top, they are
  connected in a sequence from left to right and the right most bead is connected to the
  first bead in the linker. Each bead is connected with a bond R_b = 0.47 and K_b 1250.
  There are no angles between the different head beads, but there is an angle
  [last head bead / first linker bead / first bead of first tail] Theta_a = 180 and
  K_a = 25 that helps orient the head with respect to the rest of the lipid. If an empty
  string is provided the lipid will have no head and starts with the linker beads. Spaces
  should separate the different beads; extra spaces are ignored.

  head bead types supported:
    C = NC3 = Choline      - bead Q0, charge +1
    E = NH3 = Ethanolamine - bead Qd, charge +1
    G = GL0 = Glycerol     - bead P4, charge  0
    S = CNO = Serine       - bead P5, charge  0
    P = PO4 = Phosphate    - bead Qa, charge -1
    
  Examples of lipid heads:
    "C P" -> 'NC3 PO4' - PC - PhosphatidylCholine
    "E P" -> 'NH3 PO4' - PE - PhosphatidylEthanolamine
    "G P" -> 'GLO PO4' - PG - PhosphatidylGlycerol
    "S P" -> 'CNO PO4' - PS - PhosphatidylSerine 
    "P"   -> 'PO4 ---' - PA - Phosphatidic acid
    ""    -> '--- ---' - DG - No head, Diacyl Glycerols if x2 Gly linkers are used

Linkers (-li):
  Currently only lipids with Glycerols linkers are supported. Each linker is connected
  with a bond R_b = 0.37 and K_b 1250. The number of linkers and tails provided have to
  match and each linker is connected with its corresponding tail with a bond R_b = 0.47
  and K_b = 1250. Additionally if more than one linker is provided an angle
  [last head bead / first linker bead / second linker bead] Theta_a = 120 and K_a = 25 is
  added to support the head / linker / tail orientation. 

  linker beads types supported:
    G = GLY = Glycerols    - bead Na, charge  0

  Examples of lipid heads:
    "G"     -> 'GLY --- ---' - x1 Glycerols linkers
    "G G"   -> 'GLY GLY ---' - x2 Glycerols linkers
    "G G G" -> 'GLY GLY GLY' - x3 Glycerols linkers

Tails (-ta):
  One lipid tail definition should be provided for each linker, separated with a space;
  extra spaces are ignored. Each tail can have an arbitrary number of tail beads. Tails
  are connected by bonds, first bead to the tail's linker bead and then from left to
  right (R_b = 0.47 and K_b 1250). To fix the tails orientation/dynamics an angle is
  added for each tail bead [tail bead - 1 or linker bead / tail bead / tail bead + 1]
  the Theta_a and K_a depend on the tail definition.  

  tail bead types supported:
    C = corresponding roughly to a linear combination of x4 CH2 groups (CH2-CH2-CH2-CH2).
        Represented with a C1 bead. Angle [X / C / X] Theta_a = 180 and K_a = 25.
    D = corresponding roughly to a linear combination of x2 CH2 and x2 CH groups
        (CH2-CH=CH-CH2), where the double bound is in a cis bond. Represented with a
        C3 bead, except if there is another D before or after then use a C4 bead. For
        the angles the standard [X / D / X] is a Theta_a = 120 and K_a = 45, except if
        the next bead is also D [X / D / D] then use Theta_a = 100 and K_a = 10.
    T = Same as D except with a trans double bond. The angle [X / T / X] is set with
        equilibrium bond angle Theta_a = 180 and K_a = 45. Represented with a C3 bead.
        
  Examples of tails:
    Lyso tails:
    "CCCC         " - 15-18:0 Lyso
    "CCCCC        " - 18-21:0 Lyso
    "CCCCCC       " - 21-24:0 Lyso
    Saturated tails:
    "CC     CC    " - C08-11:0 - DH 
    "CCC    CCC   " - C12-15:0 - DL or DM 
    "CCCC   CCCC  " - C15-18:0 - DP or DS
    "CCCCC  CCCCC " - C18-21:0 - DS  
    "CCCCCC CCCCCC" - C21-24:0 - DT
    Unsaturated tails:
    "CCD    CCD   " - C14:1(9c) - diMyristoleoyl 
    "CCDC   CCDC  " - C16:1(9c) - diPalmitoleoyl
    "CDDC   CDDC  " - C18:2(9c,12c) - DL or DU, di-Linoleoyl 
    "CDDD   CDDD  " - C18:3(9c,12c,15c) - diLinolenoyl
    "CCDCC  CCDCC " - C20:1(11c) - diEicosenoyl / C18:1(9c) - DO, diOleoyl 
    "DDDDC  DDDDC " - C20:4(5c,8c,11c,14c) - DA, diArachidonoyl 
    "DDDDDD DDDDDD" - C22:6(4c,7c,10c,13c,16c,19c) - diDocosahexaenoyl
    "CCCDCC CCCDCC" - C24:1(15c) - diNervonoyl
    Mixed tails:
    "CCCC   CCC   " - C14:0/16:0 - MP / C14:0/18:0 - MS
    "CCDCC  CCCC  " - C16:0/18:1(9c) - PO / C18:0/18:1(9c) - SO
    "CDDC   CCCC  " - C16:0/18:2(9c,12c) / C18:0/18:2(9c,12c)
    "DDDDC  CCCC  " - C16:0/20:4(5c,8c,11c,14c) - PA / C18:0/20:4(5c,8c,11c,14c) - SA
    "DDDDDD CCCC  " - C16:0/22:6(4c,7c,10c,13c,16c,19c) / C18:0/22:6(4c,7c,10c,13c,16c,19c)
    Trans tails:
    "CCTCC  CCTCC " - C18:1(9t) - dielaidoyl

    NOTE: the first tail is connected to GLY-1 closer to head, which is reverse order
    compared to how regular lipid names are written.

Use:
  ./lipid-martini-itp-v01.py -he 'C P' -li 'G G' -ta "CCDCC CCCC" -o POPC-lipid.itp -name POPC
"""

# Options
options = [
"""
Options:""",
("-o",       Option(str,    1,        "Martini-lipid.itp", "Output speciffic Martini lipid topology")),
("-name",    Option(str,    1,        "POPC", "Four letter lipid name")),
("-he",      Option(str,    1,        "NC3 PO4", "Lipid heads, see description")),
("-li",      Option(str,    1,        "GLY GLY", "Lipid linkers, see description")),
("-ta",      Option(str,    1,        "CCDCC CCCC", "Lipid tails, see description"))
          ]

# Define supported lipid head beads
# Lists all supported head bead types. One letter name mapped to type, atom name and charge
headMapp = {
    "C":  ['Q0', 'NC3', '1.0'], # NC3 = Choline
    "E":  ['Qd', 'NH3', '1.0'], # NH3 = Ethanolamine 
    "G":  ['P4', 'GL0', '0.0'], # GL0 = Glycerol
    "S":  ['P5', 'CNO', '0.0'], # CNO = Serine
    "P":  ['Qa', 'PO4', '-1.0'] # PO4 = Phosphate
    }

# Define possible bond lengths and forces
defBlength = '0.47'
defShortBlength = '0.37'
defBforce = '1250'

# Define possible angles and forces
defAngle1 = '100.0'
defAngle2 = '120.0'
defAngle3 = '180.0'
defAforce1 = '10.0'
defAforce2 = '25.0'
defAforce3 = '45.0'


# Get arguments    
args = sys.argv[1:]

# Print help 
if '-h' in args or '--help' in args:
    print "\n",__file__
    print desc
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
    
# Get ouput .itp file name
itpFileName  = options["-o"].value

# Get lipid description
lipidHead  = options["-he"].value
lipidLinker  = options["-li"].value
lipidTail  = options["-ta"].value
if lipidLinker==None or lipidLinker==None or lipidTail==None:
    print >>sys.stderr, "You have to provide a header, linker and tail lipid description, if one should be missing provide an empty string"
    sys.exit()
lipidName = options["-name"].value

progString = "Martini lipid itp generator version " + version + "  Args are: -o %s, -name %s, -he %s, -li %s, -ta %s" % (itpFileName, lipidName, lipidHead, lipidLinker, lipidTail)
print progString

# Open the output file
itpFile = open(itpFileName,"w")
headsArray = lipidHead.split()
linkersArray = lipidLinker.split()
linkersIndex = []
tailsArray = lipidTail.split()
if len(linkersArray)!=len(tailsArray):
    print >>sys.stderr, "A tail definition has to be provided for each linker"
    sys.exit()

# Make .itp headder
print >>itpFile, ';;;;;; Martini lipid topology auto generated using:'
print >>itpFile, '; ' + progString + '\n;\n'
print >>itpFile, '[moleculetype]'
print >>itpFile, '; molname      nrexcl'
print >>itpFile, '   ' + lipidName + '          1'

# Make lipd def
print >>itpFile, '\n[atoms]'
print >>itpFile, '; id 	type 	resnr 	residu 	atom 	cgnr 	charge'

bondsArray = []
anglesArray = []

index = 1
for cHead in headsArray:
    print >>itpFile, '%i \t%s \t%i \t%s \t%s \t%i \t%s' % (index, headMapp[cHead][0], 1, lipidName, headMapp[cHead][1], index, headMapp[cHead][2])
    if index > 1: # link head beads
        bondsArray.append([index - 1, index, defBlength, defBforce])
    index += 1
# End headsArray loop

for linkerBeadIndex in range(0, len(linkersArray)):
    cLinker = linkersArray[linkerBeadIndex]
    if cLinker != "G": # GLY
        print >>sys.stderr, "This script currently only supports GLY linkers"
        sys.exit()
    
    print >>itpFile, '%i \t%s \t%i \t%s \t%s \t%i \t%i' % (index, 'Na', 1, lipidName, "GL"+str(linkerBeadIndex+1), index, 0)
    
    linkersIndex.append(index)
    if index > 1: # There have been beads before (don't do anything if this was the first linker bead and not head)
        if len(headsArray) == (index -1):  # This is first linker bead add bond to head bead
            bondsArray.append([index - 1, index, defBlength, defBforce])
        else: # add linker linker bond
            bondsArray.append([index - 1, index, defShortBlength, defBforce])
    index += 1
# End linkersArray loop

# Orient lipid head, add angles between linkers and head + head-first linker-first tail
if len(headsArray) > 0:
    # Add angle between the last head bead and first x2 linkers
    if len(linkersArray) > 1:
        anglesArray.append([linkersIndex[0] - 1, linkersIndex[0], linkersIndex[1], defAngle2, defAforce2])
    # Add angle between the last head bead, the first linker and the first tail
    anglesArray.append([linkersIndex[0] - 1, linkersIndex[0], linkersIndex[-1]+1, defAngle3, defAforce2])


tailIndex = 0
indexToLetter = "A B C D E F G H I J K L M N".split()
for cTail in tailsArray:   
    # Add bond from current tail to assosiated linker
    bondsArray.append([linkersIndex[tailIndex], index, defBlength, defBforce])

    for cTailBeadIndex in range(0, len(cTail)):
        cTailBead = cTail[cTailBeadIndex]
        if cTailBead=='C':
            bType = 'C1'
        elif cTailBead=='T':
            bType = 'C3'
        elif cTailBead=='D':
            if (cTailBeadIndex > 0 and cTail[cTailBeadIndex - 1] == 'D') or (len(cTail) > (cTailBeadIndex + 1) and cTail[cTailBeadIndex + 1] == 'D'):
                # Another D before or after
                bType = 'C4'
            else:
                bType = 'C3'
        else:
            print >>sys.stderr, "Tail definition \"%s\" not recognized" % (cTailBead)
            sys.exit()
        
        atomName = cTailBead + str(cTailBeadIndex+1) + indexToLetter[tailIndex]       
        print >>itpFile, '%i \t%s \t%i \t%s \t%s \t%i \t%i' % (index, bType, 1, lipidName, atomName, index, 0)

        # Add bond between tail beads
        if cTailBeadIndex > 0:
            bondsArray.append([index-1, index, defBlength, defBforce])

        # Add angles to support the tails (only add angles when considdering the middle bead in the angle)
        if (cTailBeadIndex + 1) < len(cTail): # Else we are looking at the last bead (which can't have an angle)
            # Set angle constreins (default regular 'C' angle)
            cdefAngle = defAngle3
            cdefAforce = defAforce2
            if cTail[cTailBeadIndex]=='D': 
                if cTail[cTailBeadIndex + 1]=='D':  # has another D right after
                    cdefAngle = defAngle1
                    cdefAforce = defAforce1                    
                else: # does not have another D after
                    cdefAngle = defAngle2
                    cdefAforce = defAforce3
            elif cTail[cTailBeadIndex]=='T':  # Trans double bond
                cdefAngle = defAngle3
                cdefAforce = defAforce3

            if cTailBeadIndex == 0: # top angle connecting to the tail linker [linker / head tail / head+1 tail]
                anglesArray.append([linkersIndex[tailIndex], index, index + 1, cdefAngle, cdefAforce])
            else: # regular angle connecting to the tail linker [current-1 tail / current tail / current+1 tail bead]
                anglesArray.append([index - 1, index, index + 1, cdefAngle, cdefAforce])             
        # end angle stuff
         
        index += 1
    tailIndex += 1
# End tailsArray loop    

# Write lipid bonds
print >>itpFile, '\n[bonds]'
print >>itpFile, '; i j 	funct 	length 	force.c.'
for cBond in bondsArray:
    print >>itpFile, '  %i %i\t1 \t%s \t%s' % (cBond[0], cBond[1], cBond[2], cBond[3])

# Write lipid angles
print >>itpFile, '\n[angles]'
print >>itpFile, '; i j k 	funct 	angle 	force.c.'
for cAngle in anglesArray:
    print >>itpFile, '  %i %i %i \t2 \t%s \t%s' % (cAngle[0], cAngle[1], cAngle[2], cAngle[3], cAngle[4])

itpFile.close()
# End lipid-martini-itp



