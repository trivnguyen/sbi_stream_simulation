
from ml_collections import config_dict

def get_config():
    config = config_dict.ConfigDict()

    # working directory
    config.workdir = './output'
    config.binary_dir = './binaries'
    config.potential_dir = './pot'
    config.save_observables_only = False
    config.overwrite = True

    # simulation configuration
    config.num_samples = 2  # set to a small number for testing

    # prior configuration
    config.priors = priors = config_dict.ConfigDict()
    priors.M_sat = dict(dist='log_uniform', a=1e-5, b=1e-1) # 1e10 Msun
    priors.rs_sat = dict(dist='log_uniform', a=1e-2, b=3)  # kpc  # match NFW profile
    priors.vz = dict(dist='normal', loc=0, scale=100)  # km/s
    priors.vphi = dict(dist='normal', loc=0, scale=100)  # km/s
    priors.r = dict(dist='uniform', loc=0, scale=3)  # 0 - few kpc  # 3 = high
    priors.phi = dict(dist='uniform', loc=0, scale=360.)  # TODO: periodic parameters
    # TODO: change phi_a to use 90th percentile stream length
    priors.phi_a = dict(dist='uniform', loc=-40, scale=80)  # impact point w.r.t prog
    priors.t_a = dict(dist='uniform', loc=0, scale=5)  # kpc/km/s  # difficult one # 5 = integration time
    priors.tmax = dict(dist='delta', value=4)  # integration time, leave at 4

    return config
