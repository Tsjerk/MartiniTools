\#!/usr/bin/env python

import sys

version = "1.1"
tag     = "20120826-TAW"


# Martini Force Field Definition Generating Script - BMW version
# Version 2012-02-19 TAW


# If a file with nonbonded parameters (atom types and interactions)
# is given on the command line, the output will be a merged forcefield
aa = len(sys.argv) > 1 and open(sys.argv[1]).readlines()


beads = { # Name: Mass
    "P5":   72,  "P4":   72,  "P3":   72,  "P2":   72,  "P1":   72, # Polar
    "Nda":  72,  "Nd":   72,  "Na":   72,  "N0":   72,              # Intermediate
    "C5":   72,  "C4":   72,  "C3":   72,  "C2":   72,  "C1":   72, # Apolar
    "Qda":  72,  "Qd":   72,  "Qa":   72,  "Q0":   72,              # Charged
    "SP5":  45,  "SP4":  45,  "SP3":  45,  "SP2":  45,  "SP1":  45, # Ring type, polar
    "SNda": 45,  "SNd":  45,  "SNa":  45,  "SN0":  45,              # Ring type, intermediate
    "SC5":  45,  "SC4":  45,  "SC3":  45,  "SC2":  45,  "SC1":  45, # Ring type, apolar
    "SQda": 45,  "SQd":  45,  "SQa":  45,  "SQ0":  45,              # Ring type, charged
    "AC1":  72,  "AC2":  72,                                        # Amino acid side chains (Q-C interactions)
    "RQd":  72,  "AQa":  72,                                        # BMW charged types
    "BMC":  72,  "BMQ":  72,                                        # BMW water beads    
}

martini_v2_1 ="""
; MARTINI FORCEFIELD V2.1 - FINAL VERSION WITH ADAPTATIONS FOR BMW WATER MODEL
;
; SJ MARRINK, 16-06-2007 (last modified: 26-1-2011)
;
; please cite:
;
; L. Monticelli, S. Kandasamy, X. Periole, R. Larson, D.P. Tieleman, S.J. Marrink.
; The MARTINI coarse grained force field: extension to proteins.
; J. Chem. Th. Comp., 4:819-834, 2008. 
;
; S.J. Marrink, H.J. Risselada, S. Yefimov, D.P. Tieleman, A.H. de Vries.
; The MARTINI forcefield: coarse grained model for biomolecular simulations.
; JPC-B, 111:7812-7824, 2007.
;
; and (if using lipid topologies):
;
; S.J. Marrink, A.H. de Vries, A.E. Mark.
; Coarse grained model for semi-quantitative lipid simulations.
; JPC-B, 108:750-760, 2004.


[ defaults ]
1 1

[ atomtypes ]

; Currently eighteen main particle types are defined, 
; divided into four main categories:
; (P, polar; N, intermediate; C, apolar; Q, charged)
; each of which has a number of sublevels (0,a,d, or ad) 
; subtype 0 has no hydrogen bond forming capacities,
; subtype d has hydrogen donor capacities, 
; subtype a has hydrogen acceptor capacities, 
; and subtype da has both donor and acceptor capacities 
; or (1,2,3,4,5) where subtype 5 is more polar than 1.

; Two main classes of particles are furthermore distinguished, namely
; STANDARD particles which are mapped using a 4:1 mapping scheme, and
; RING particles which are used for ring compounds mapped 2-3:1.
; A special BIG particle type is defined in addition to prevent freezing of CG water.
; Two AMINO acid particle types are used to avoid Q-C clashes inside proteins.

; For reasons of computational efficiency, all particle masses are set to 72 amu,
; except for ring types which are set to 45 amu.
; For realistic dynamics, the particle masses should be adapted. 
; This might require a reduction of the integration timestep, however.

; name mass charge ptype c6 c12
"""


# Epsilon: 5.60(A), 5.00(B), 4.50(C), 4.00(D), 3.50(E), 3.10(F), 2.70(G), 2.30(H), 2.00(I)
# Sigma:   0.00(0), 0.43(1), 0.47(2), 0.57(3), 0.62(4)
# Scaling: 1.00(a), 0.95(b), 0.90(c), 0.75(d), 0.71(e)

epsilon = {
    "0":0.00, # no interaction
    "A":5.60, # supra attractive
    "B":5.00, # attractive
    "C":4.50, # almost attractive
    "D":4.00, # semi attractive
    "E":3.50, # intermediate
    "F":3.10, # almost intermediate
    "G":2.70, # semi repulsive
    "H":2.30, # almost repulsive
    "I":2.00, # repulsive
    "Z":0.25, # with sigma=1 yields c6=12=1
    }

sigma   = {
    "0":0.00, # no interaction
    "1":0.43, # ring bead type
    "2":0.47, # default bead type 
    "3":0.57, # supra attractive bead type
    "4":0.62, # super repulsive bead type
    "5":1.00, # unity: c6/12 parity
    }

scale   = {
    "a":1.00, # no scaling
    "b":0.95, # 95%: interaction with polarizable water
    "c":0.90, # 90%: S* - C1 interactions with BMW water
    "d":0.75, # 75%: ring bead types (S*) and uncharged interactions with BMW for epsilon >=4.5
    "e":0.71, # 71%: uncharged interactions with BMW for epsilon<4.5
    }

# BEAD TYPES

# Default bead types; 4:1 mapping, 72 AMU
plain = { 
    "P5":   72,  "P4":   72,  "P3":   72,  "P2":   72,  "P1":   72, # Polar
    "Nda":  72,  "Nd":   72,  "Na":   72,  "N0":   72,              # Intermediate
    "C5":   72,  "C4":   72,  "C3":   72,  "C2":   72,  "C1":   72, # Apolar
    "Qda":  72,  "Qd":   72,  "Qa":   72,  "Q0":   72,              # Charged
    "AC1":  72,  "AC2":  72,                                        # Amino acid specific (Q-C interactions)
    }

# Ring bead types; mapping 2:1 or 3:1, 45 AMU
small = {
    "SP5":  45,  "SP4":  45,  "SP3":  45,  "SP2":  45,  "SP1":  45, # Ring type, polar
    "SNda": 45,  "SNd":  45,  "SNa":  45,  "SN0":  45,              # Ring type, intermediate
    "SC5":  45,  "SC4":  45,  "SC3":  45,  "SC2":  45,  "SC1":  45, # Ring type, apolar
    "SQda": 45,  "SQd":  45,  "SQa":  45,  "SQ0":  45,              # Ring type, charged
    }

# Virtual sites, plain type; mapping 4:1, no mass
vsite = {
    "vP5":   0,  "vP4":   0,  "vP3":   0,  "vP2":   0,  "vP1":   0, # Polar
    "vNda":  0,  "vNd":   0,  "vNa":   0,  "vN0":   0,              # Intermediate
    "vC5":   0,  "vC4":   0,  "vC3":   0,  "vC2":   0,  "vC1":   0, # Apolar
    "vQda":  0,  "vQd":   0,  "vQa":   0,  "vQ0":   0,              # Charged
    "vAC1":  0,  "vAC2":  0,                                        # Amino acid specific (Q-C interactions)
    }

# Virtual sites, small type; mapping 2:1 or 3:1, no mass
svste = {
    "vSP5":  0,  "vSP4":  0,  "vSP3":  0,  "vSP2":  0,  "vSP1":  0, # Ring type, polar
    "vSNda": 0,  "vSNd":  0,  "vSNa":  0,  "vSN0":  0,              # Ring type, intermediate
    "vSC5":  0,  "vSC4":  0,  "vSC3":  0,  "vSC2":  0,  "vSC1":  0, # Ring type, apolar
    "vSQda": 0,  "vSQd":  0,  "vSQa":  0,  "vSQ0":  0,              # Ring type, charged
    }

# Other: Special purpose types
other = {
    "BP4":  72, # Big water particle (antifreeze). Not used with BMW
    "D":     0, # Dummy particle type
    "BMC":  72, # BMW water model central bead
    "BMQ":   0, # BMW water model satellite
    "RQd":  72, # BMW guanidinium group (Arginine)
    "AQa":  72, # BMW negative bead
    }

classes = ("plain","small","vsite","svste","other")

# Dummy atom types. These will be given a repulsion "DUMMY_REPEL"
# with all atoms from the atomistic forcefield, if provided
dummy = ["D","BMQ"]

# Gather all atom types
all,mass = zip(*[i for j in classes for i in eval(j).items()])
virtual  = vsite.keys) + svste.keys()

rla   = range(len(all))
cmb   = [ (all[i],all[j]) for i in rla for j in rla[i:] ]

table_plain = """
       Qda   Qd   Qa   Q0   P5   P4   P3   P2   P1  Nda   Nd   Na   N0   C5   C4   C3   C2   C1
  Qda   A    A    A    C    A    A    A    B    B    A    A    A    E    D    E    E    E    E 
   Qd   A    B    A    C    A    A    A    B    B    A    D    A    E    D    E    E    E    E 
   Qa   A    A    B    C    A    A    A    B    B    A    A    D    E    D    E    E    E    E 
   Q0   C    C    C    E    B    A    B    C    D    D    D    D    E    D    E    E    E    E 
   P5   A    A    A    B    A    A    A    A    A    A    A    A    E    F    G    G    H    I 
   P4   A    A    A    A    A    B    B    C    C    D    A    A    A    F    G    G    H    I 
   P3   A    A    A    B    A    B    B    C    C    C    C    C    E    E    F    F    G    H 
   P2   B    B    B    C    A    C    C    C    C    C    C    C    D    E    E    F    G    H 
   P1   B    B    B    D    A    C    C    C    C    C    C    C    D    E    E    E    F    G 
  Nda   A    A    A    D    A    A    C    C    C    C    C    C    E    E    F    G    G    G 
   Nd   A    D    A    D    A    A    C    C    C    C    D    C    E    E    F    G    G    G 
   Na   A    A    D    D    A    A    C    C    C    C    C    D    E    E    F    G    G    G 
   N0   E    E    E    E    E    E    E    D    D    E    E    E    E    E    E    E    F    G 
   C5   D    D    D    D    F    F    E    E    E    E    E    E    E    E    E    E    F    F 
   C4   E    E    E    E    G    G    F    E    E    F    F    F    E    E    E    E    F    F 
   C3   E    E    E    E    G    G    F    F    E    G    G    G    E    E    E    E    E    E 
   C2   E    E    E    E    H    H    G    G    F    G    G    G    F    F    F    E    E    E 
   C1   E    E    E    E    I    I    H    H    G    G    G    G    G    F    F    E    E    E
"""

table_small = """
      SQda  SQd  SQa  SQ0  SP5  SP4  SP3  SP2  SP1 SNda  SNd  SNa  SN0  SC5  SC4  SC3  SC2  SC1
 SQda  A1d  A1d  A1d  C1d  A1d  A1d  A1d  B1d  B1d  A1d  A1d  A1d  E1d  D1d  E1d  E1d  E1d  E1d
  SQd  A1d  B1d  A1d  C1d  A1d  A1d  A1d  B1d  B1d  A1d  D1d  A1d  E1d  D1d  E1d  E1d  E1d  E1d
  SQa  A1d  A1d  B1d  C1d  A1d  A1d  A1d  B1d  B1d  A1d  A1d  D1d  E1d  D1d  E1d  E1d  E1d  E1d
  SQ0  C1d  C1d  C1d  E1d  B1d  A1d  B1d  C1d  D1d  D1d  D1d  D1d  E1d  D1d  E1d  E1d  E1d  E1d
  SP5  A1d  A1d  A1d  B1d  A1d  A1d  A1d  A1d  A1d  B1d  B1d  B1d  E1d  F1d  G1d  G1d  H1d  I1d
  SP4  A1d  A1d  A1d  A1d  A1d  B1d  B1d  C1d  C1d  D1d  D1d  D1d  A1d  F1d  G1d  G1d  H1d  I1d
  SP3  A1d  A1d  A1d  B1d  A1d  B1d  B1d  C1d  C1d  C1d  C1d  C1d  E1d  E1d  F1d  F1d  G1d  H1d
  SP2  B1d  B1d  B1d  C1d  A1d  C1d  C1d  C1d  C1d  C1d  C1d  C1d  D1d  E1d  E1d  F1d  G1d  H1d
  SP1  B1d  B1d  B1d  D1d  A1d  C1d  C1d  C1d  C1d  C1d  C1d  C1d  D1d  E1d  E1d  E1d  F1d  G1d
 SNda  A1d  A1d  A1d  D1d  B1d  D1d  C1d  C1d  C1d  C1d  C1d  C1d  E1d  E1d  F1d  G1d  G1d  G1d
  SNd  A1d  D1d  A1d  D1d  B1d  D1d  C1d  C1d  C1d  C1d  D1d  C1d  E1d  E1d  F1d  G1d  G1d  G1d
  SNa  A1d  A1d  D1d  D1d  B1d  D1d  C1d  C1d  C1d  C1d  C1d  D1d  E1d  E1d  F1d  G1d  G1d  G1d
  SN0  E1d  E1d  E1d  E1d  E1d  E1d  E1d  D1d  D1d  E1d  E1d  E1d  E1d  E1d  E1d  E1d  F1d  G1d
  SC5  D1d  D1d  D1d  D1d  F1d  F1d  E1d  E1d  E1d  E1d  E1d  E1d  E1d  E1d  E1d  E1d  F1d  F1d
  SC4  E1d  E1d  E1d  E1d  G1d  G1d  F1d  E1d  E1d  F1d  F1d  F1d  E1d  E1d  E1d  E1d  F1d  F1d
  SC3  E1d  E1d  E1d  E1d  G1d  G1d  F1d  F1d  E1d  G1d  G1d  G1d  E1d  E1d  E1d  E1d  E1d  E1d
  SC2  E1d  E1d  E1d  E1d  H1d  H1d  G1d  G1d  F1d  G1d  G1d  G1d  F1d  F1d  F1d  E1d  E1d  E1d
  SC1  E1d  E1d  E1d  E1d  I1d  I1d  H1d  H1d  G1d  G1d  G1d  G1d  G1d  F1d  F1d  E1d  E1d  E1d
"""

table_small_plain = """
       Qda   Qd   Qa   Q0   P5   P4   P3   P2   P1  Nda   Nd   Na   N0   C5   C4   C3   C2   C1
 SQda   A    A    A    C    A    A    A    B    B    A    A    A    E    D    E    E    E    Ec
  SQd   A    B    A    C    A    A    A    B    B    A    D    A    E    D    E    E    E    Ec
  SQa   A    A    B    C    A    A    A    B    B    A    A    D    E    D    E    E    E    Ec
  SQ0   C    C    C    E    B    A    B    C    D    D    D    D    E    D    E    E    E    Ec
  SP5   A    A    A    B    A    A    A    A    A    A    A    A    E    F    G    G    H    Ic
  SP4   A    A    A    A    A    B    B    C    C    D    A    A    A    F    G    G    H    Ic
  SP3   A    A    A    B    A    B    B    C    C    C    C    C    E    E    F    F    G    Hc
  SP2   B    B    B    C    A    C    C    C    C    C    C    C    D    E    E    F    G    Hc
  SP1   B    B    B    D    A    C    C    C    C    C    C    C    D    E    E    E    F    Gc
 SNda   A    A    A    D    A    A    C    C    C    C    C    C    E    E    F    G    G    Gc
  SNd   A    D    A    D    A    A    C    C    C    C    D    C    E    E    F    G    G    Gc
  SNa   A    A    D    D    A    A    C    C    C    C    C    D    E    E    F    G    G    Gc
  SN0   E    E    E    E    E    E    E    D    D    E    E    E    E    E    E    E    F    Gc
  SC5   D    D    D    D    F    F    E    E    E    E    E    E    E    E    E    E    F    Fc
  SC4   E    E    E    E    G    G    F    E    E    F    F    F    E    E    E    E    F    Fc
  SC3   E    E    E    E    G    G    F    F    E    G    G    G    E    E    E    E    E    Ec
  SC2   E    E    E    E    H    H    G    G    F    G    G    G    F    F    F    E    E    Ec
  SC1   E    E    E    E    I    I    H    H    G    G    G    G    G    F    F    E    E    Ec
"""

table_vsite_plain = """
       Qda   Qd   Qa   Q0   P5   P4   P3   P2   P1  Nda   Nd   Na   N0   C5   C4   C3   C2   C1
 vQda   A    A    A    C    A    A    A    B    B    A    A    A    E    D    E    E    E    E 
  vQd   A    B    A    C    A    A    A    B    B    A    D    A    E    D    E    E    E    E 
  vQa   A    A    B    C    A    A    A    B    B    A    A    D    E    D    E    E    E    E 
  vQ0   C    C    C    E    B    A    B    C    D    D    D    D    E    D    E    E    E    E 
  vP5   A    A    A    B    A    A    A    A    A    A    A    A    E    F    G    G    H    I 
  vP4   A    A    A    A    A    B    B    C    C    D    A    A    A    F    G    G    H    I 
  vP3   A    A    A    B    A    B    B    C    C    C    C    C    E    E    F    F    G    H 
  vP2   B    B    B    C    A    C    C    C    C    C    C    C    D    E    E    F    G    H 
  vP1   B    B    B    D    A    C    C    C    C    C    C    C    D    E    E    E    F    G 
 vNda   A    A    A    D    A    A    C    C    C    C    C    C    E    E    F    G    G    G 
  vNd   A    D    A    D    A    A    C    C    C    C    D    C    E    E    F    G    G    G 
  vNa   A    A    D    D    A    A    C    C    C    C    C    D    E    E    F    G    G    G 
  vN0   E    E    E    E    E    E    E    D    D    E    E    E    E    E    E    E    F    G 
  vC5   D    D    D    D    F    F    E    E    E    E    E    E    E    E    E    E    F    F 
  vC4   E    E    E    E    G    G    F    E    E    F    F    F    E    E    E    E    F    F 
  vC3   E    E    E    E    G    G    F    F    E    G    G    G    E    E    E    E    E    E 
  vC2   E    E    E    E    H    H    G    G    F    G    G    G    F    F    F    E    E    E 
  vC1   E    E    E    E    I    I    H    H    G    G    G    G    G    F    F    E    E    E
"""

table_small_vsite = """
      vQda  vQd  vQa  vQ0  vP5  vP4  vP3  vP2  vP1 vNda  vNd  vNa  vN0  vC5  vC4  vC3  vC2  vC1
 SQda   A    A    A    C    A    A    A    B    B    A    A    A    E    D    E    E    E    Ec
  SQd   A    B    A    C    A    A    A    B    B    A    D    A    E    D    E    E    E    Ec
  SQa   A    A    B    C    A    A    A    B    B    A    A    D    E    D    E    E    E    Ec
  SQ0   C    C    C    E    B    A    B    C    D    D    D    D    E    D    E    E    E    Ec
  SP5   A    A    A    B    A    A    A    A    A    A    A    A    E    F    G    G    H    Ic
  SP4   A    A    A    A    A    B    B    C    C    D    A    A    A    F    G    G    H    Ic
  SP3   A    A    A    B    A    B    B    C    C    C    C    C    E    E    F    F    G    Hc
  SP2   B    B    B    C    A    C    C    C    C    C    C    C    D    E    E    F    G    Hc
  SP1   B    B    B    D    A    C    C    C    C    C    C    C    D    E    E    E    F    Gc
 SNda   A    A    A    D    A    A    C    C    C    C    C    C    E    E    F    G    G    Gc
  SNd   A    D    A    D    A    A    C    C    C    C    D    C    E    E    F    G    G    Gc
  SNa   A    A    D    D    A    A    C    C    C    C    C    D    E    E    F    G    G    Gc
  SN0   E    E    E    E    E    E    E    D    D    E    E    E    E    E    E    E    F    Gc
  SC5   D    D    D    D    F    F    E    E    E    E    E    E    E    E    E    E    F    Fc
  SC4   E    E    E    E    G    G    F    E    E    F    F    F    E    E    E    E    F    Fc
  SC3   E    E    E    E    G    G    F    F    E    G    G    G    E    E    E    E    E    Ec
  SC2   E    E    E    E    H    H    G    G    F    G    G    G    F    F    F    E    E    Ec
  SC1   E    E    E    E    I    I    H    H    G    G    G    G    G    F    F    E    E    Ec
"""

table_svste_plain = """
       Qda   Qd   Qa   Q0   P5   P4   P3   P2   P1  Nda   Nd   Na   N0   C5   C4   C3   C2   C1
vSQda   A    A    A    C    A    A    A    B    B    A    A    A    E    D    E    E    E    Ec
 vSQd   A    B    A    C    A    A    A    B    B    A    D    A    E    D    E    E    E    Ec
 vSQa   A    A    B    C    A    A    A    B    B    A    A    D    E    D    E    E    E    Ec
 vSQ0   C    C    C    E    B    A    B    C    D    D    D    D    E    D    E    E    E    Ec
 vSP5   A    A    A    B    A    A    A    A    A    A    A    A    E    F    G    G    H    Ic
 vSP4   A    A    A    A    A    B    B    C    C    D    A    A    A    F    G    G    H    Ic
 vSP3   A    A    A    B    A    B    B    C    C    C    C    C    E    E    F    F    G    Hc
 vSP2   B    B    B    C    A    C    C    C    C    C    C    C    D    E    E    F    G    Hc
 vSP1   B    B    B    D    A    C    C    C    C    C    C    C    D    E    E    E    F    Gc
vSNda   A    A    A    D    A    A    C    C    C    C    C    C    E    E    F    G    G    Gc
 vSNd   A    D    A    D    A    A    C    C    C    C    D    C    E    E    F    G    G    Gc
 vSNa   A    A    D    D    A    A    C    C    C    C    C    D    E    E    F    G    G    Gc
 vSN0   E    E    E    E    E    E    E    D    D    E    E    E    E    E    E    E    F    Gc
 vSC5   D    D    D    D    F    F    E    E    E    E    E    E    E    E    E    E    F    Fc
 vSC4   E    E    E    E    G    G    F    E    E    F    F    F    E    E    E    E    F    Fc
 vSC3   E    E    E    E    G    G    F    F    E    G    G    G    E    E    E    E    E    Ec
 vSC2   E    E    E    E    H    H    G    G    F    G    G    G    F    F    F    E    E    Ec
 vSC1   E    E    E    E    I    I    H    H    G    G    G    G    G    F    F    E    E    Ec
"""

table_svste_small = """
      SQda  SQd  SQa  SQ0  SP5  SP4  SP3  SP2  SP1 SNda  SNd  SNa  SN0  SC5  SC4  SC3  SC2  SC1
vSQda  A1d  A1d  A1d  C1d  A1d  A1d  A1d  B1d  B1d  A1d  A1d  A1d  E1d  D1d  E1d  E1d  E1d  E1d
 vSQd  A1d  B1d  A1d  C1d  A1d  A1d  A1d  B1d  B1d  A1d  D1d  A1d  E1d  D1d  E1d  E1d  E1d  E1d
 vSQa  A1d  A1d  B1d  C1d  A1d  A1d  A1d  B1d  B1d  A1d  A1d  D1d  E1d  D1d  E1d  E1d  E1d  E1d
 vSQ0  C1d  C1d  C1d  E1d  B1d  A1d  B1d  C1d  D1d  D1d  D1d  D1d  E1d  D1d  E1d  E1d  E1d  E1d
 vSP5  A1d  A1d  A1d  B1d  A1d  A1d  A1d  A1d  A1d  B1d  B1d  B1d  E1d  F1d  G1d  G1d  H1d  I1d
 vSP4  A1d  A1d  A1d  A1d  A1d  B1d  B1d  C1d  C1d  D1d  D1d  D1d  A1d  F1d  G1d  G1d  H1d  I1d
 vSP3  A1d  A1d  A1d  B1d  A1d  B1d  B1d  C1d  C1d  C1d  C1d  C1d  E1d  E1d  F1d  F1d  G1d  H1d
 vSP2  B1d  B1d  B1d  C1d  A1d  C1d  C1d  C1d  C1d  C1d  C1d  C1d  D1d  E1d  E1d  F1d  G1d  H1d
 vSP1  B1d  B1d  B1d  D1d  A1d  C1d  C1d  C1d  C1d  C1d  C1d  C1d  D1d  E1d  E1d  E1d  F1d  G1d
vSNda  A1d  A1d  A1d  D1d  B1d  D1d  C1d  C1d  C1d  C1d  C1d  C1d  E1d  E1d  F1d  G1d  G1d  G1d
 vSNd  A1d  D1d  A1d  D1d  B1d  D1d  C1d  C1d  C1d  C1d  D1d  C1d  E1d  E1d  F1d  G1d  G1d  G1d
 vSNa  A1d  A1d  D1d  D1d  B1d  D1d  C1d  C1d  C1d  C1d  C1d  D1d  E1d  E1d  F1d  G1d  G1d  G1d
 vSN0  E1d  E1d  E1d  E1d  E1d  E1d  E1d  D1d  D1d  E1d  E1d  E1d  E1d  E1d  E1d  E1d  F1d  G1d
 vSC5  D1d  D1d  D1d  D1d  F1d  F1d  E1d  E1d  E1d  E1d  E1d  E1d  E1d  E1d  E1d  E1d  F1d  F1d
 vSC4  E1d  E1d  E1d  E1d  G1d  G1d  F1d  E1d  E1d  F1d  F1d  F1d  E1d  E1d  E1d  E1d  F1d  F1d
 vSC3  E1d  E1d  E1d  E1d  G1d  G1d  F1d  F1d  E1d  G1d  G1d  G1d  E1d  E1d  E1d  E1d  E1d  E1d
 vSC2  E1d  E1d  E1d  E1d  H1d  H1d  G1d  G1d  F1d  G1d  G1d  G1d  F1d  F1d  F1d  E1d  E1d  E1d
 vSC1  E1d  E1d  E1d  E1d  I1d  I1d  H1d  H1d  G1d  G1d  G1d  G1d  G1d  F1d  F1d  E1d  E1d  E1d
"""

table_other = """
        D  BMC  BMQ  RQd  AQa
    D   0    0    0    0    0 
  BMC   0    Z5   0    E    B
  BMQ   0    0    0    0    0
  RQd   0    E    0    B    A
  AQa   0    B    0    A    B
"""

table_other_plain = """
       Qda   Qd   Qa   Q0   P5   P4   P3   P2   P1  Nda   Nd   Na   N0   C5   C4   C3   C2   C1
    D   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0 
  BMC   B    B    B    B    Ad   Bd   Bd   Cd   Cd   De   De   De   Ee   Fe   Ge   Ge   He   Ie
  BMQ   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0 
  RQd   A    B    A    C    A    A    A    B    B    A    D    A    E    D    E    E    E    E
  AQa   A    A    B    C    A    A    A    B    B    A    A    D    E    D    A    A    A    A
"""

table_other_small = """
      SQda  SQd  SQa  SQ0  SP5  SP4  SP3  SP2  SP1 SNda  SNd  SNa  SN0  SC5  SC4  SC3  SC2  SC1
    D   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0 
  BMC   B    B    B    B    Ad   Bd   Bd   Cd   Cd   De   De   De   Ee   Fe   Ge   Ge   He   Ie
  BMQ   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0 
  RQd   A    B    A    C    A    A    A    B    B    A    D    A    E    D    E    E    E    E
  AQa   A    A    B    C    A    A    A    B    B    A    A    D    E    D    A    A    A    A
"""

table_other_vsite = """
      vQda  vQd  vQa  vQ0  vP5  vP4  vP3  vP2  vP1 vNda  vNd  vNa  vN0  vC5  vC4  vC3  vC2  vC1
    D   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0 
  BMC   B    B    B    B    Ad   Bd   Bd   Cd   Cd   De   De   De   Ee   Fe   Ge   Ge   He   Ie
  BMQ   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0 
  RQd   A    B    A    C    A    A    A    B    B    A    D    A    E    D    E    E    E    E
  AQa   A    A    B    C    A    A    A    B    B    A    A    D    E    D    A    A    A    A
"""

table_other_svste = """
      vSQda vSQd vSQa vSQ0 vSP5 vSP4 vSP3 vSP2 vSP1 vSNda vSNd vSNa vSN0 vSC5 vSC4 vSC3 vSC2 vSC1
    D   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0 
  BMC   B    B    B    B    Ad   Bd   Bd   Cd   Cd   De   De   De   Ee   Fe   Ge   Ge   He   Ie
  BMQ   0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0 
  RQd   A    B    A    C    A    A    A    B    B    A    D    A    E    D    E    E    E    E
  AQa   A    A    B    C    A    A    A    B    B    A    A    D    E    D    A    A    A    A
"""


def table2pairs(x):
    labels2 = x.pop(0)
    return [((i[0],j),k) for i in x for j,k in zip(labels2,i[1:])]

def sigeps2c(eps=None,scl="a",sig="2"):
    """
    Convert a string encoding for epsilon, sigma and scaling to C6 and C12 parameters
    """
    if eps == None: 
        return None,None
    if scl in sigma.keys(): 
        sig, scl = scl, "2"
    return 4*epsilon[eps]*scale[scl]*sigma[sig]**6, 4*epsilon[eps]*scale[scl]*sigma[sig]**12
    
tables = []
for i in classes:
    tables.append(globals().get("table_%s"%i))
    for j in classes:
        tables.append(globals().get("table_%s_%s"%(i,j)))
tables = [[ j.split() for j in i.split("\n") if j not in ("","\n") ] for i in tables if i]

pairs  = dict([j for i in tables for j in table2pairs(i)])


## -- Prepare the output
atomtypes      = []
nonbond_params = []
pairtypes      = []


# Read in atomistic stuff
if len(sys.argv) > 1:
    aa         = open(sys.argv[1]).readlines()
    key        = ""

    atomtypes.append("; Atomistic definitions\n")
    nonbond_params.append("; Atomistic definitions\n")
    pairtypes.append("; Atomistic definitions\n")

    for line in aa:
        s = line.strip()
        if s and s[0] == "[":
            key=line            
        elif "nonbond_params" in key:
            nonbond_params.append(line)
        elif "atomtypes" in key:
            atomtypes.append(line)
        elif "pairtypes" in key:
            pairtypes.append(line)
    atomtypes.append("; End of atomistic definitions\n")
    nonbond_params.append("; End of atomistic definitions\n")
    pairtypes.append("; End of atomistic definitions\n")

# Note which atomtypes we have for defining interactions with dummy particles             
aa_atomtypes = [i.split()[0] for i in atomtypes if i.strip() and not i.strip()[0] == ";"]

# Add coarsegrained stuff
typestr="%5s    0 %10.3f      0.000     %1s   0.0           0.0\n"
atomtypes.extend([typestr%(tp,ms,tp in virtual and 'V' or 'A') for tp,ms in zip(all,mass)])

for i,j in cmb:
    c6,c12 = sigeps2c(*pairs.get((i,j),pairs.get((j,i),"")))
    if c6:
        nonbond_params.append(" %7s  %7s %2d  %e %e\n"%(i,j,1,c6,c12))

atomtypes.append("; End of coarsegrained definitions\n")
if nonbond_params:
    nonbond_params.append("; End of coarsegrained definitions\n")
if pairtypes:
    pairtypes.append("; End of coarsegrained definitions\n")

for i in dummy:
    if i in tp:
        for j in aa_atomtypes:
            nonbond_params.append(" %7s  %7s  %2d 0.0 DUMMY_REPEL\n"%(i,j,1))

# Print stuff:

print "; This file was created automagically by", sys.argv[0] 
print "; (c)2012 Tsjerk A Wassenaar, University of Groningen"
print ";"

if len(sys.argv) > 1:
    print "; This file contains a merged forcefield, combining %s with MARTINI" % sys.argv[1]
    print ";"
    print "#define DUMMY_REPEL 1e-7"
    print ";"

print martini_v2_1

print "".join(atomtypes),      "\n"
print "[ nonbond_params ]\n", "".join(nonbond_params), "\n"
if pairtypes:
    print "[ pairtypes ]\n",      "".join(pairtypes),      "\n"
print

for i in sys.argv[2:]: print open(i).read()
    


    
