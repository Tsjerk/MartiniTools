import versions 

from tprio import *

# Separate modules:
import adress       # AdResS (Adaptive Resolution Simulation)
import efield       # Electric Field
import fep          # Free Energy Perturbation
import imd          # Interactive Molecular Dynamics
import pull         # Pull Groups
import rotate       # Enforced rotation 
import elphys       # Electric Physiology 
import qmmm         # QM/MM 
import wall         # Walls

        
class Opts(ListWithNames):
    def __init__(self,tpr):        
        self.position = tpr.tell()

        self.extend([
            ("ngtc",              Integer(tpr)),                                # Number of temperature groups
            ("nh_chain_length",   Integer(tpr) if tpr.version >= 69 else None), # Nose-Hoover chain length
            ("ngacc",             Integer(tpr)),                                # Number of acceleration groups
            ("ngfrz",             Integer(tpr)),                                # Number of freeze groups
            ("ngener",            Integer(tpr)),                                # Number of energy groups
        ])

        self.extend([ 
            ("nrdf",              Tuple(tpr, self.ngtc,      Integer if tpr.version < 13 else Real )),  # Numbers of degrees of freedom 
            ("ref_t",             Tuple(tpr, self.ngtc,      Real )),                                   # Reference temperatures
            ("tau_t",             Tuple(tpr, self.ngtc,      Real )),                                   # Coupling times
            ("nFreeze",           Tuple(tpr, self.ngfrz,     (Integer, Integer, Integer) )),            # Freeze dimensions (x,y,z) per group
            ("acc",               Tuple(tpr, self.ngacc,     (Real, Real, Real) )),                     # Accelerations (x,y,z) per group
            ("egp_flags",         Tuple(tpr, self.ngener**2, Integer ) if tpr.version >= 12 else None), # 
            ("annealing",         Tuple(tpr, self.ngtc,      Integer ) if tpr.version >  26 else None), # Annealing type per group
            ("annealing_npoints", Tuple(tpr, self.ngtc,      Integer ) if tpr.version >  26 else ()),   # Number of points per group
            ])

        self.extend([ 
            ("anneal_time",       Tuple(tpr, self.annealing_npoints, Real)), # Times between annealing points
            ("anneal_temp",       Tuple(tpr, self.annealing_npoints, Real)), # Temperatures at annealing points
            ])

        # print "Temperature coupling:", zip(self.ref_t, self.tau_t)
        # print "Freeze groups ({}): {} {}".format(self.nFreeze.position, self.ngfrz, self.nFreeze)
        # print "Acceleration groups ({}): {} {}".format(self.acc.position, self.ngacc, self.acc)
        # print "Anneal times ({}):        {}".format(self.anneal_time.position, self.anneal_time)
        # print "Anneal temperatures ({}): {}".format(self.anneal_temp.position, self.anneal_temp)
            
        
        
class InputRecord(ListWithNames):
    def __init__(self,tpr):

        self.extend([
            # The MD/EM/... integrator; index to enum
            ("integrator",           Integer(tpr)),

            # Number of steps to run
            ("nsteps",               Integer(tpr)    if tpr.version < 62         else Long(tpr)),  

            # Starting step of run (part)
            ("init_step",            Integer(tpr)    if tpr.version < 62         else Long(tpr)),  

            # Number of simulation part 
            ("simulation_part",      Integer(tpr)    if tpr.version >= 58        else None),       

            # Frequency of calculating energies
            ("nstcalcenergy",        Integer(tpr)    if tpr.version >= 67        else None),       

            # Type of periodic boundary conditions; index to enum
            ("pbc",                  Integer(tpr)    if tpr.version <  53        else None),

            # Molecules wrapping over PBC no/yes
            ("periodic_molecules",   Integer(tpr)    if (45 <= tpr.version < 53) else None),

            # Cutoff scheme; index to enum
            ("cutoff_scheme",        Integer(tpr)    if tpr.version >= 80        else None),

            # Neighbour search type; index to enum
            ("ns_type",              Integer(tpr)),
            
            # Frequency of neighbour list update
            ("nstlist",              Integer(tpr)),

            # Number of cells per rlong
            ("ndelta",               Integer(tpr)),

            
            ("dummy",                Integer(tpr)    if tpr.version <  41        else None),
            ("dummy",                Integer(tpr)    if tpr.version <  41        else None),
            
            # Test particle insertion radius
            ("rtpi",                 Real(tpr)       if tpr.version >= 45        else None),

            ("nstcomm",              Integer(tpr)),
            ("comm_mode",            Integer(tpr)    if tpr.version >  34        else None),
            ("nstcheckpoint",        Integer(tpr)    if tpr.version >  25 and tpr.version < versions.RemoveObsoleteParameters1 else None), 

            # Frequency of steepest descents steps when doing CG
            ("nstcgsteep",           Integer(tpr)),

            # NBFG-S minimization
            ("nbfgscorr",            Integer(tpr)    if tpr.version >= 30        else None),

            # Frequency of writing to log
            ("nstlog",               Integer(tpr)),

            # Frequency of writing coordinates to TRR
            ("nstxout",              Integer(tpr)),
            
            # Frequency of writing velocities to TRR
            ("nstvout",              Integer(tpr)),

            # Frequency of writing forces to TRR
            ("nstfout",              Integer(tpr)),

            # Frequency of writing energies to EDR
            ("nstenergy",            Integer(tpr)),

            # Frequency of writing coordinates to XTC
            ("nstxtcout",            Integer(tpr)),

            # Starting time
            ("tinit",                Real(tpr)       if tpr.version <  59        else Double(tpr)), 

            # Time step
            ("dt",                   Real(tpr)       if tpr.version <  59        else Double(tpr)), 

            # Precision of XTC file 
            ("xtc_precision",        Real(tpr)),                                                

            ("dummy",                Integer(tpr)    if tpr.version <  19        else None),
            ("dummy",                Integer(tpr)    if tpr.version <  19        else None),
            ("dummy",                Integer(tpr)    if tpr.version <  18        else None),

            # Allowed energy drift due to the Verlet buffer in kJ/mol/ps per atom,
            # a value of -1 means: use rlist
            ("verlet_buffer_drift",  Real(tpr)       if tpr.version >= 80        else None),

            # Neighbour list cut-off
            ("rlist",                Real(tpr)),

            # Long-range cut-off for switched potentials
            ("rlistlong",            Real(tpr)       if tpr.version >= 67        else None),

            # Frequency of long-range cut-off calculation
            ("nstcalclr",            Integer(tpr)    if tpr.version >= 82 and tpr.version != 90 else None),      

            # Treatment of electrostatic interactions; index to enum
            ("coulombtype",          Integer(tpr)),

            # Modifier for electrostatic interactions; index to enum
            ("coulomb_modifier",     Integer(tpr)    if tpr.version >= 81        else None),

            # Radius to start switching coulomb interactions
            ("rcoulomb_switch",      Real(tpr)),
            
            # Cut-off radius for coulomb interactions
            ("rcoulomb",             Real(tpr)),

            # Treatment of vanderwaals interactions; index to enum
            ("vdwtype",              Integer(tpr)),

            # Modifier for vanderwaals interactions; index to enum
            ("vdw_modifier",         Integer(tpr)    if tpr.version >= 81        else None),

            # Radius to start switching vanderwaals interactions
            ("rvdw_switch",          Real(tpr)),

            # Cut-off radius for vanderwaals interactions
            ("rvdw",                 Real(tpr)),

            # Apply long range dispersion corrections for Energy and Pressure
            ("eDispCorr",            Integer(tpr)),

            # Dielectric constant of medium
            ("epsilon_r",            Real(tpr)),

            # Dielectric constant of reaction field
            ("epsilon_rf",           Real(tpr)       if tpr.version >= 37        else None),

            # Extension of the potential lookup tables beyond the cut-off
            ("table_extension",      Real(tpr)       if tpr.version >= 29        else None),
            
            # Algorithm for calculating Born radii            
            ("gb_algorithm",         Integer(tpr)    if tpr.version >= 25        else None),

            # Frequency of calculating the Born radii inside rlist
            ("nstgbradii",           Integer(tpr)    if tpr.version >= 25        else None),

            # Cutoff for Born radii calculation# the contribution from atoms
            # between rlist and rgbradii is updated every nstlist steps
            ("rgbradii",             Real(tpr)       if tpr.version >= 25        else None),

            # Salt concentration in M for Generalized Born models
            ("gb_saltconc",          Real(tpr)       if tpr.version >= 25        else None),

            # Implicit solvent algorithm; index to enum.
            ("implicit_solvent",     Integer(tpr)    if tpr.version >= 25        else None),

            # Dielectric coefficient of the implicit solvent
            ("gb_epsilon_solvent",   Real(tpr)       if tpr.version >= 55        else None),

            # Scaling factors used in the OBC GB model. Default values are OBC(II)
            ("gb_obc_alpha",         Real(tpr)       if tpr.version >= 55        else None),
            ("gb_obc_beta",          Real(tpr)       if tpr.version >= 55        else None),
            ("gb_obc_gamma",         Real(tpr)       if tpr.version >= 55        else None),
            ("gb_dielectric_offset", Real(tpr)       if tpr.version >= 60        else None),

            # Algorithm for non-polar part of GBSA; index to enum
            ("sa_algorithm",         Integer(tpr)    if tpr.version >= 60        else None),

            # Surface tension (kJ/mol/nm^2) for the SA (nonpolar surface) part of GBSA
            # The value -1 will set default value for Still/HCT/OBC GB-models.
            ("sa_surface_tension",   Real(tpr)       if tpr.version >= 55        else None),

            # Spacing for the PME/PPPM FFT grid
            ("fourierspacing",       Real(tpr)       if tpr.version >= 80        else None),

            # FFT grid size, when a value is 0 fourierspacing will be used
            ("fourier_nx",           Integer(tpr) ),
            ("fourier_ny",           Integer(tpr) ),
            ("fourier_nz",           Integer(tpr) ),

            # EWALD/PME/PPPM parameters
            ("pme_order",            Integer(tpr) ),
            ("ewald_rtol",           Real(tpr) ),
            ("ewald_rtol_lj",        Real(tpr)       if tpr.version >= 93        else None),
            # Index to enum
            ("ewald_geometry",       Integer(tpr)    if tpr.version >= 24        else None), 

            ("dummy",                Integer(tpr)    if tpr.version == 17        else None), 
            ("epsilon_surface",      Real(tpr)       if tpr.version  > 17        else None),

            # Boolean for optimization of FFT settings
            ("optimize_fft",         Integer(tpr)    if tpr.version < versions.RemoveObsoleteParameters1 else None),

            ("ljpmeCombRule",        Integer(tpr)    if tpr.version >= 93        else None),

            # Boolean for continuation of run (no constraints in first step)
            ("continuation",         Integer(tpr) ),
            
            # Temperature coupling algorithm; index to enum
            ("tcoupl",               Integer(tpr) ),

            # ...
            ("print_nose_hoover_chain_variables",       Integer(tpr)    if tpr.version >= 79        else None),  

            # Frequency of applying temperature coupling correction
            ("nsttcouple",           Integer(tpr)    if tpr.version >= 71        else None),  

            ("dummy",                Integer(tpr)    if tpr.version <= 15        else None), 

            ("epct",                 Integer(tpr)    if tpr.version <= 17        else None),
            ("dummy",                Integer(tpr)    if tpr.version <= 15        else None),

            # Pressure coupling method; index to enum
            ("pcoupl",               Integer(tpr)    if tpr.version >  17        else None),

            # Type of pressure coupling (Isotropic, Semiisotropci, Anisotropic); index to enum
            ("pcoupltype",           Integer(tpr)    if tpr.version >  17        else None),

            # Frequency of pressure coupling
            ("nstpcouple",           Integer(tpr)    if tpr.version >= 71        else None),

            # Relaxation time for pressure coupling
            ("tau_p",                Real(tpr) ),

            # Reference pressure
            ("ref_p",                RealVector(tpr) if tpr.version <= 15        else RealMatrix(tpr)),

            # Compressibility of the system
            ("compressibility",      RealVector(tpr) if tpr.version <= 15        else RealMatrix(tpr)),

            # Scaling of reference coordinates with pressure coupling; index to enum
            ("refcoord_scaling",     Integer(tpr)    if tpr.version >= 47        else None), 


            
            ("posres_com",           RealVector(tpr) if tpr.version >= 47        else None),
            ("posres_comB",          RealVector(tpr) if tpr.version >= 47        else None),
            ("andersen_seed",        Integer(tpr)    if (25 < tpr.version < 79)  else None),
            ("bSimAnn",              Integer(tpr)    if tpr.version <  26        else None),
            ("zerotemptime",         Real(tpr)       if tpr.version <  26        else None),
            ("dummy",                Real(tpr)       if tpr.version <  37        else None),

            # Relative tolerance of SHAKE 
            ("shake_tol",            Real(tpr)),

            
            ("fudgeQQ",              Real(tpr)       if tpr.version <  54        else None),
            ("efep",                 fep.FEP(tpr)),
            ])

        self.extend([
            ("simtemp",              tpr.version >= 79 and fep.SimTemp(tpr,self.efep.n_lambda)),
            ("expanded",             tpr.version >= 79 and fep.Expanded(tpr,self.efep)),

            # Distance restraints
            ("disre",                Integer(tpr)    if tpr.version >= 57        else None ),
            ("disre_weighting",      Integer(tpr)),
            ("bDisreMixed",          Integer(tpr)),
            ("disre_fc",             Real(tpr)),
            ("disre_tau",            Real(tpr)),
            ("nstdisreout",          Integer(tpr)),

            # Orientation restraints
            ("orire_fc",             Real(tpr)       if tpr.version >= 22        else None),                           #
            ("orire_tau",            Real(tpr)       if tpr.version >= 22        else None),                           #
            ("nstorireout",          Integer(tpr)    if tpr.version >= 22        else None),                           #

            # Dihedral restraint force constant
            ("dihre_fc",             Real(tpr)       if 26 <= tpr.version < 79   else None),                           #

            ("dummy",                Real(tpr)       if 26 <= tpr.version < 56   else None),                           #
            ("dummy",                Integer(tpr)    if 26 <= tpr.version < 56   else None),                           #

            # Initial stepsize for EM
            ("emstep",               Real(tpr)),

            # Force tolerance for EM
            ("emtol",                Integer(tpr)),
            
            ("bShakeSOR",            Integer(tpr)    if tpr.version >= 22        else None),

            # Max number of iterations in relax-shells
            ("niter",                Integer(tpr)    if tpr.version >= 11        else None),

            # Step size (ps^2) for minimization of flexible constraints
            ("fcstep",               Real(tpr)       if tpr.version >= 21        else None),

            # Constraint algorithm; index to enum
            ("eConstrAlg",           Integer(tpr)),

            # Projection order of LINCS 
            ("lincs_order",          Integer(tpr)),                                                                    #

            # LINCS warns if angle changes more than this
            ("lincs_warnangle",      Real(tpr)),                                                                       #

            ("dummy",                Integer(tpr)    if tpr.version <= 14        else None),                           #

            # Number of iterations in final step of LINCS
            ("lincs_iter",           Integer(tpr)    if tpr.version >= 26        else None),                           #

            # ...
            ("bd_temp",              Real(tpr)       if tpr.version <  33        else None),                           #

            # Brownian Dynamics friction
            ("bd_fric",              Real(tpr)),                                                                       

            # Langevin Dynamics random seed
            ("ld_seed",              Long(tpr)       if tpr.version >= versions.Use64BitRandomSeed else Integer(tpr)), 

            # Deformation matrix (Non-equilibrium MD)
            ("deform",               RealMatrix(tpr) if tpr.version >= 33        else None),

            # Cosine shaped acceleration (Non-equilibrium MD)
            ("cos_acceleration",     Real(tpr)       if tpr.version >= 14        else None),

            # User definable parameters
            ("userint1",             Integer(tpr)),
            ("userint2",             Integer(tpr)),
            ("userint3",             Integer(tpr)),
            ("userint4",             Integer(tpr)),
            ("userreal1",            Float(tpr, tpr.precision)),
            ("userreal2",            Float(tpr, tpr.precision)),
            ("userreal3",            Float(tpr, tpr.precision)),
            ("userreal4",            Float(tpr, tpr.precision)),
        ])

        self.extend([
            ("adress",               adress.AdResS(tpr)),
            ("pull",                 pull.Pull(tpr)),
            ("rot",                  rotate.Rotate(tpr)),
            ("imd",                  imd.InteractiveMD(tpr)),
            ("opts",                 Opts(tpr)),
            ("wall",                 wall.Wall(tpr)),
            ("efield",               efield.Efield(tpr)),
            ("swapcoords",           elphys.SwapCoords(tpr)),
            ("qmmm",                 qmmm.QMMM(tpr)),
        ])

        # print self.lincs_warnangle

        # print "nsteps", self.nsteps
        # print "nstlist", self.nstlist
        # print "nstcomm", self.nstcomm
        # print "comm_mode", self.comm_mode
        # print "nstcheckpoint", self.nstcheckpoint
        # print "nstcgsteep", self.nstcgsteep
        # print "nbfgscorr", self.nbfgscorr
        # print "nstlog", self.nstlog
        # print "nstxout", self.nstxout
        # print "nstvout", self.nstvout
        # print "nstfout", self.nstfout
        # print "nstenergy", self.nstenergy
        # print "nstxtcout", self.nstxtcout
        # print "init_t", self.t_init, self.t_init.position
        # print "delta_t", self.dt
        # print "xtcprec", self.xtcprec
        # print "rvdw", self.rvdw
        # print "fourier_spacing", self.fourier_spacing
        # print "bContinuation", self.bContinuation
        # print "nsttcouple", self.nsttcouple
        # print "nstpcouple", self.nstpcouple
        # print "Pressure coupling time", self.tau_p
        # print "Reference pressure:\n", self.ref_p
        # print "Compressibility:\n", self.compressibility
        # print "SHAKE tolerance:", self.shake_tol
        # print "EFEP:", hasattr(self.efep,"empty") and "NO" or self.efep
        # print "SIMTEMP:", self.simtemp
        # print "disre_fc:", self.disre_fc
        # print "disre_tau:", self.disre_tau
        # print "NSTDISREOUT:", self.nstdisreout
        # print "NSTORIREOUT:", self.nstdisreout
        # print "LINCS warning angle:", self.lincs_warnangle
        # print "Cosine acceleration:", self.cos_acceleration
        # print "User INT  1:", self.userint1
        # print "User INT  2:", self.userint2
        # print "User INT  3:", self.userint3
        # print "User INT  4:", self.userint4
        # print "User REAL 1:", self.userreal1
        # print "User REAL 2:", self.userreal2
        # print "User REAL 3:", self.userreal3
        # print "User REAL 4:", self.userreal4

        # "Th-Th-Th-Th-Th-... That's all, folks."
        
        
