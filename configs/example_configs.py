
from ml_collections import config_dict

def get_config():
    config = config_dict.ConfigDict()

    # working directory
    config.workdir = './output'
    config.binary_dir = './binaries'
    config.potential_dir = './pot'
    config.save_observables_only = False

    # simulation configuration
    config.num_samples = 2  # set to a small number for testing

    # prior configuration
    config.priors = priors = config_dict.ConfigDict()
    priors.M_sat = dict(dist='log_uniform', a=1e-5, b=1)
    priors.vz = dict(dist='normal', loc=-50, scale=50)
    priors.r = dict(dist='delta', value=0.2)
    priors.phi = dict(dist='delta', value=250.)
    priors.vphi = dict(dist='delta', value=35.)
    priors.t_a = dict(dist='delta', value=0.2)
    priors.phi_a = dict(dist='delta', value=-4)
    priors.tmax = dict(dist='delta', value=4)

    return config
