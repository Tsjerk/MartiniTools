#
#	File 'md-prod-out.mdp' was generated
#	By user: tsjerk (501)
#	On host: Tsjerk.local
#	At date: Mon Dec 29 14:36:13 2014
#

# VARIOUS PREPROCESSING OPTIONS
# Preprocessor information: use cpp syntax.
# e.g.: -I/home/joe/doe -I/home/mary/roe
StringIgnore include                  = 

# e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)
StringIgnore: define                   = 

# RUN CONTROL PARAMETERS
Enum: integrator               = md

# Start time and timestep in ps
Real if version < 59 else Double:  tinit                    = 0
Real if version < 59 else Double:  dt                       = 0.002
Integer if version < 62 else Long: nsteps                   = 2500000

# For exact run continuation or redoing part of a run
Integer if version < 62 else Long: init_step                = 0

# Part index is updated automatically on checkpointing (keeps files separate)
Integer: simulation_part          = 1

# mode for center of mass motion removal
Enum (None|Linear|Angular|RTC?): comm_mode                = Linear

# number of steps for center of mass motion removal
Integer: nstcomm                  = 100

# group(s) for center of mass motion removal
Groupnames[6]: comm_grps                = System

# LANGEVIN DYNAMICS OPTIONS
# Friction coefficient (amu/ps) and random seed
Real: bd_fric                  = 0
Long if version >= versions.Use64BitRandomSeed else Integer: ld_seed                  = 1993

# ENERGY MINIMIZATION OPTIONS
# Force tolerance and initial step-size
Integer: emtol                    = 10
Real:    emstep                   = 0.01
# Max number of iterations in relax-shells
Integer: niter                    = 20
# Step size (ps^2) for minimization of flexible constraints
Real:    fcstep                   = 0

# Frequency of steepest descents steps when doing CG
Integer: nstcgsteep               = 1000

Integer: nbfgscorr                = 10


# TEST PARTICLE INSERTION OPTIONS
Real: rtpi                     = 0.05


# OUTPUT CONTROL OPTIONS
# Output frequency for coords (x), velocities (v) and forces (f)
Integer: nstxout                  = 500
Integer: nstvout                  = 500
Integer: nstfout                  = 0

# Output frequency for energies to log file and energy file
Integer: nstlog                   = 25000
Integer: nstcalcenergy            = 100
Integer: nstenergy                = 500

# Output frequency and precision for .xtc file
Integer: nstxtcout                = 25000
Real:    xtc_precision            = 1000

# This selects the subset of atoms for the .xtc file. You can
# select multiple groups. By default all atoms will be written.
Groupnames[7]: xtc_grps                 = 

# Selection of energy groups
Groupnames[1]: energygrps               = Solute Solvent


# NEIGHBORSEARCHING PARAMETERS
# cut-off scheme (group: using charge groups, Verlet: particle based cut-offs)
Enum: cutoff_scheme            = Group

# nblist update frequency
Integer: nstlist                  = 5

# ns algorithm (simple or grid)
Enum: ns_type                  = Grid

# Periodic boundary conditions: xyz, no, xy
Enum: pbc                      = xyz

Boolean: periodic_molecules       = no

# Allowed energy drift due to the Verlet buffer in kJ/mol/ps per atom,
# a value of -1 means: use rlist
Real: verlet_buffer_drift      = 0.005

# nblist cut-off        
Real: rlist                    = 1.2

# long-range cut-off for switched potentials
Real: rlistlong                = -1

Integer: nstcalclr                = -1


# OPTIONS FOR ELECTROSTATICS AND VDW
# Method for doing electrostatics
Enum: coulombtype              = PME-switch
Enum: coulomb_modifier         = Potential-shift-Verlet
Real: rcoulomb_switch          = 0
Real: rcoulomb                 = 1.1

# Relative dielectric constant for the medium and the reaction field
Real: epsilon_r                = 1
Real: epsilon_rf               = 0

# Method for doing Van der Waals
Enum: vdwtype                  = shift
Enum: vdw_modifier             = Potential-shift-Verlet

# cut-off lengths       
Real: rvdw_switch              = 0.9
Real: rvdw                     = 1.0

# Apply long range dispersion corrections for Energy and Pressure
Enum: DispCorr                 = EnerPres

# Extension of the potential lookup tables beyond the cut-off
Real: table_extension          = 1

# Separate tables between energy group pairs
energygrp_table          = 

# Spacing for the PME/PPPM FFT grid
Real: fourierspacing           = 0.125

# FFT grid size, when a value is 0 fourierspacing will be used
Integer: fourier_nx               = 0
Integer: fourier_ny               = 0
Integer: fourier_nz               = 0

# EWALD/PME/PPPM parameters
Integer: pme_order                = 4
Real:    ewald_rtol               = 1e-05
Enum:    ewald_geometry           = 3d
Real:    epsilon_surface          = 0
Boolean: optimize_fft             = yes


# IMPLICIT SOLVENT ALGORITHM
Enum:    implicit_solvent         = No


# GENERALIZED BORN ELECTROSTATICS
# Algorithm for calculating Born radii
Enum: gb_algorithm             = Still

# Frequency of calculating the Born radii inside rlist
Integer: nstgbradii               = 1

# Cutoff for Born radii calculation# the contribution from atoms
# between rlist and rgbradii is updated every nstlist steps
Real: rgbradii                 = 1

# Dielectric coefficient of the implicit solvent
Real: gb_epsilon_solvent       = 80

# Salt concentration in M for Generalized Born models
Real: gb_saltconc              = 0

# Scaling factors used in the OBC GB model. Default values are OBC(II)
Real: gb_obc_alpha             = 1
Real: gb_obc_beta              = 0.8
Real: gb_obc_gamma             = 4.85
Real: gb_dielectric_offset     = 0.009

Enum: sa_algorithm             = Ace-approximation

# Surface tension (kJ/mol/nm^2) for the SA (nonpolar surface) part of GBSA
# The value -1 will set default value for Still/HCT/OBC GB-models.
Real: sa_surface_tension       = -1


# OPTIONS FOR WEAK COUPLING ALGORITHMS
# Temperature coupling  
Enum:    tcoupl                   = v-rescale
Integer: nsttcouple               = 5
Opts->Integer: nh_chain_length          = 10
Boolean: print_nose_hoover_chain_variables = no

# Groups to couple separately
Groupnames[0]: tc_grps                  = Solute Solvent

# Time constant (ps) and reference temperature (K)
Opts->Tuple(Real) tau_t                    =  0.1 0.1
Opts->Tuple(Real) ref_t                    = 1000 1000

# pressure coupling     
Enum:    pcoupl                   = no
Enum:    pcoupltype               = Isotropic
Integer: nstpcouple               = 5

# Time constant (ps), compressibility (1/bar) and reference P (bar)
Real:    tau_p                    = 1.0
RealMatrix: compressibility          = 4.5e-5
RealMatrix: ref_p                    = 1
# Scaling of reference coordinates, No, All or COM
Enum: refcoord_scaling         = No


# OPTIONS FOR QMMM calculations
Boolean: QMMM                     = no
# Groups treated Quantum Mechanically
Groupnames[9]: QMMM_grps                = 
# QM method             
Tuple(Enum): QMmethod                 = 
# QMMM scheme           
Tuple(Enum): QMMMscheme               = normal
# QM basisset           
Tuple(Enum): QMbasis                  = 
# QM charge             
Tuple(Integer): QMcharge                 = 
# QM multiplicity       
Tuple(Integer): QMmult                   = 
# Surface Hopping       
Tuple(Boolean): SH                       = 
# CAS space options     
Tuple(Integer): CASorbitals              = 
Tuple(Integer): CASelectrons             = 
Tuple(Real): SAon                     = 
Tuple(Real): SAoff                    = 
Tuple(Integer): SAsteps                  = 
# Scale factor for MM charges
Real: MMChargeScaleFactor      = 1
# Optimization of QM subsystem
Tuple(Boolean): bOPT                     = 
Tuple(Boolean): bTS                      = 

# SIMULATED ANNEALING  
# Type of annealing for each temperature group (no/single/periodic)
Tuple(Enum): annealing                = 
# Number of time points to use for specifying annealing in each group
Tuple(Integer): annealing_npoints        = 
# List of times at the annealing points for each group
Tuple(Group(Real)) annealing_time           = 
# Temp. at each annealing point, for each group.
Tuple(Group(Real)) annealing_temp           = 

# GENERATE VELOCITIES FOR STARTUP RUN
gen-vel                  = no
gen-temp                 = 300
gen-seed                 = 173529

# OPTIONS FOR BONDS    
constraints              = all-bonds
# Type of constraint algorithm
Enum: constraint_algorithm     = Lincs
# Do not constrain the start configuration
Boolean: continuation             = no

# Use successive overrelaxation to reduce the number of shake iterations
Shake_SOR                = no

# Relative tolerance of shake
Real: shake_tol                = 0.0001

# Highest order in the expansion of the constraint coupling matrix
Integer: lincs_order              = 4

# Number of iterations in the final step of LINCS. 1 is fine for
# normal simulations, but use 2 to conserve energy in NVE runs.
# For energy minimization with constraints it should be 4 to 8.
Integer: lincs_iter               = 2

# Lincs will write a warning to the stderr if in one step a bond
# rotates over more degrees than
Real: lincs_warnangle          = 30

# Convert harmonic bonds to morse potentials
morse                    = no

# ENERGY GROUP EXCLUSIONS
# Pairs of energy groups for which all non-bonded interactions are excluded
energygrp_excl           = 

# WALLS                
# Number of walls, type, atom types, densities and box-z scale factor for Ewald
nwall                    = 0
wall-type                = 9-3
wall-r-linpot            = -1
wall-atomtype            = 
wall-density             = 
wall-ewald-zfac          = 3

# COM PULLING          
# Pull type: no, umbrella, constraint or constant-force
pull                     = no

# ENFORCED ROTATION    
# Enforced rotation: No or Yes
rotation                 = no

# NMR refinement stuff 
# Distance restraints type: No, Simple or Ensemble
Enum: disre                    = No
# Force weighting of pairs in one distance restraint: Conservative or Equal
Enum: disre_weighting          = Conservative
# Use sqrt of the time averaged times the instantaneous violation
Boolean: disre_mixed              = no
Real: disre_fc                 = 1000
Real: disre_tau                = 0
# Output frequency for pair distances to energy file
Integer: nstdisreout              = 100

# Orientation restraints: No or Yes
orire                    = no
# Orientation restraints force constant and tau for time averaging
Real: orire_fc                 = 0
Real: orire_tau                = 0
orire_fitgrp             = 
# Output frequency for trace(SD) and S to energy file
Integer: nstorireout              = 100

# Free energy variables
Boolean: free_energy              = no

couple_moltype           = 
couple_lambda0           = vdw-q
couple_lambda1           = vdw-q
couple_intramol          = no

Double:  init_lambda              = -1
Integer: init_lambda_state        = -1
Double:  delta_lambda             = 0

Integer: nstdhdl                  = 50

Tuple(Double): fep_lambdas              = 
Tuple(Double): mass_lambdas             = 
Tuple(Double): coul_lambdas             = 
Tuple(Double): vdw_lambdas              = 
Tuple(Double): bonded_lambdas           = 
Tuple(Double): restraint_lambdas        = 
Tuple(Double): temperature_lambdas      = 

Integer: calc_lambda_neighbors    = 1

init_lambda_weights      = 

Boolean: dhdl_print_energy        = no
Real: sc_alpha                 = 0
Integer: sc_power                 = 1
Real: sc_r_power               = 6
Real: sc_sigma                 = 0.3
Boolean: sc_coul                  = no
Boolean: separate_dhdl_file       = yes
Boolean: dhdl_derivatives         = yes
Integer: dh_hist_size             = 0
Double: dh_hist_spacing          = 0.1

# Non-equilibrium MD stuff
acc-grps                 = 
accelerate               = 
freezegrps               = 
freezedim                = 
Real: cos_acceleration         = 0
RealMatrix: deform                   = 

# simulated tempering variables
SimTemp: simulated_tempering      = no
Enum: simulated_tempering_scaling = geometric
Real: sim_temp_low             = 300
Real: sim_temp_high            = 300

# Electric fields      
# Format is number of terms (int) and for all terms an amplitude (real)
# and a phase angle (real)
Efield: E_x                      = 
Efield: E_xt                     = 
Efield: E_y                      = 
Efield: E_yt                     = 
Efield: E_z                      = 
Efield: E_zt                     = 

# AdResS parameters    
adress                   = no

# User defined thingies
Groupnames[4]: user1_grps               = 
Groupnames[5]: user2_grps               = 
Integer: userint1                 = 0
Integer: userint2                 = 0
Integer: userint3                 = 0
Integer: userint4                 = 0
Real: userreal1                = 0
Real: userreal2                = 0
Real: userreal3                = 0
Real: userreal4                = 0
