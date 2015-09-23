import versions

from tprio import *

class SimTemp(ListWithNames):
    def __init__(self,tpr,n_lambda):
        self.append(("simulated_tempering", Integer(tpr)))

        # Method; index to enum
        self.extend([
            ("simulated_tempering_scaling", Integer(tpr) if self.simulated_tempering else None),
            ("sim_temp_high", Real(tpr) if self.simulated_tempering else None),
            ("sim_temp_low",  Real(tpr) if self.simulated_tempering else None),
            ("temperatures",  Tuple(tpr, n_lambda, Real) if self.simulated_tempering else None),
        ])
        
class Expanded(versions.Test):
    def __init__(self,tpr,fep):
        if not Integer(tpr):
            return
        self.init_lambda_weights = fep.n_lambda and [ Float(tpr, tpr.precision) for i in range(fep.n_lambda) ]
        self.bInit_weights       = fep.n_lambda and Integer(tpr)
        self.nstexpanded         = Integer(tpr)
        self.elmcmove            = Integer(tpr)
        self.elamstats           = Integer(tpr)
        self.lmc_repeats         = Integer(tpr)
        self.gibbsdeltalam       = Integer(tpr)
        self.lmc_forced_nstart   = Integer(tpr)
        self.lmc_seed            = Integer(tpr)
        self.mc_temp             = Float(tpr, tpr.precision)
        self.bSymmetrizedTMatrix = Integer(tpr)
        self.nstTij              = Integer(tpr)
        self.minvarmin           = Integer(tpr)
        self.c_range             = Integer(tpr)
        self.wl_scale            = Float(tpr, tpr.precision)
        self.wl_ratio            = Float(tpr, tpr.precision)
        self.init_wl_delta       = Float(tpr, tpr.precision)
        self.bWLoneovert         = Integer(tpr)
        self.elmceq              = Integer(tpr)
        self.equil_steps         = Integer(tpr)
        self.equil_samples       = Integer(tpr)
        self.equil_n_at_lam      = Integer(tpr)
        self.equil_wl_delta      = Float(tpr, tpr.precision)
        self.equil_ratio         = Float(tpr, tpr.precision)


class FEP(ListWithNames):
    types = ["fep_lambdas", "mass_lambdas", "coul_lambdas", "vdw_lambdas", 
             "bonded_lambdas", "restraint_lambdas", "temperature_lambdas"]

    def __init__(self,tpr):
        # Free energy calculations, yes/no
        self.extend([
            ("free_energy", Integer(tpr)),
            ("init_lambda_state",  Integer(tpr) if tpr.version >= 79 else None),
            ("init_lambda",        Double(tpr) if tpr.version >= 59 else Real(tpr)),        
            ("delta_lambda",       Double(tpr) if tpr.version >= 59 else Real(tpr)), 
            ("n_lambda",           Integer(tpr) if tpr.version >= 64 else 0),
        ])

        self.all_lambda     = []
        if tpr.version >= 79:
            for typ in self.types:
                self.append((typ, Tuple(tpr, self.n_lambda, Double)))
                self.separate_dvdl = self.n_lambda and Tuple(tpr, len(self.types), Integer) 
            if self.n_lambda == 0 and self.init_lambda >= 0:
                self.separate_dvdl = [True]
        elif tpr.version >= 64:
            self.all_lambda.append(Tuple(tpr, self.n_lambda, Double))


        # Soft-core potentials
        self.extend([
            ("sc_alpha",              Real(tpr)    if tpr.version >= 13 else None),       
            ("sc_power",              Integer(tpr) if tpr.version >= 38 else None),       
            ("sc_r_power",            Real(tpr)    if tpr.version >= 79 else 6.0),
            ("sc_sigma",              Real(tpr)    if tpr.version >= 15 else 0.3),
            ("sc_coul",               Integer(tpr) if tpr.version >= 79 else None),       
            
            # Output frequency of dH/dL
            ("nstdhdl",               Integer(tpr) if tpr.version >= 64 else None),       
            
            #
            ("separate_dhdl_file",    Integer(tpr) if tpr.version >= 73 else None),       
            ("dhdl_derivatives",      Integer(tpr) if tpr.version >= 73 else None),       
            ("dh_hist_size",          Integer(tpr) if tpr.version >= 71 else None),       
            ("dh_hist_spacing",       Double(tpr)  if tpr.version >= 71 else None),       
            ("dhdl_print_energy",     Integer(tpr) if tpr.version >= 79 else None),     
            
            ("calc_lambda_neighbors", Integer(tpr) if tpr.version >= 83 else None),
        ])
