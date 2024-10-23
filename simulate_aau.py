
import os
import shutil
import sys
import h5py
import time
from typing import Optional, Tuple, Union

sys.path.append('/global/u2/t/tvnguyen/sbi_stream_simulation/sims')

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
from tqdm import tqdm
from absl import flags, logging
from ml_collections import config_flags, config_dict

import initial_stream
import stream_impact
import subhalo_orbit
from rotation_matrix import obs_from_pos6d

logging.set_verbosity(logging.INFO)

def read_posteriors(filename, num_samples=None, quantiles=None, seed=None):
    """ Read in the posteriors and return the parameters """
    # read data from posterior.hdf5
    with h5py.File(filename, 'r') as f:
        print(f.keys())
        M_sat = f['M_sat'][:]
        rs_sat = f['rs_sat'][:]
        vz = f['vz'][:]
        vphi = f['vphi'][:]
        r = f['r'][:]
        phi = f['phi'][:]
        t_a = f['t_a'][:]
        phi_a = f['phi_a'][:]
        tmax = f['tmax'][:]

    # subsample to num_samples
    num_data_tot = len(M_sat)
    if num_samples is not None:
        rng = np.random.RandomState(seed)
        subset = rng.choice(num_data_tot, num_samples, replace=False)
        M_sat = M_sat[subset]
        rs_sat = rs_sat[subset]
        vz = vz[subset]
        vphi = vphi[subset]
        r = r[subset]
        phi = phi[subset]
        t_a = t_a[subset]
        phi_a = phi_a[subset]
        tmax = tmax[subset]
    num_data = len(M_sat)

    if quantiles is not None:
        M_sat = np.nanquantile(M_sat, quantiles)
        rs_sat = np.nanquantile(rs_sat, quantiles)
        vz = np.nanquantile(vz, quantiles)
        vphi = np.nanquantile(vphi, quantiles)
        r = np.nanquantile(r, quantiles)
        phi = np.nanquantile(phi, quantiles)
        t_a = np.nanquantile(t_a, quantiles)
        phi_a = np.nanquantile(phi_a, quantiles)
        tmax = np.nanquantile(tmax, quantiles)

    params_dict = dict(
        r=r, phi=phi, vphi=vphi, vz=vz, M_sat=M_sat, tmax=tmax, t_a=t_a,
        phi_a=phi_a, rs_sat=rs_sat
    )

    return params_dict, num_data


def run_single_sim(
    r, phi, vphi, vz, M_sat, tmax, t_a, phi_a, rs_sat, pid):
    """ Simulate a stream with a subhalo impact.

    Parameters
    ----------
    r : float
        Impact parameter in kpc (distance from stream to sat)
    phi : float
        Angle around stream in dec
    vphi : float
        Velocity around stream in km/s
    vz : float
        Velocity along stream in km/s
    M_sat : float
        Mass of subhalo in 1e10 Msun
    tmax : float
        How long stream disrupts in Gyr
    t_a : float
        Time since interaction in Gyr
    phi_a : float
        Interaction point along stream in deg, phi=0 is progenitor location
        (try -20 to 10)
    rs_sat : float
        Scale radius of subhalo in kpc, can be adjusted along with M_sat
        using equation 15 in erkal et al. 2015
    pid : int
        Index included in saved filenames

    Returns
    -------
    """
    # TODO: this hack needs to be streamlined better
    DEFAULT_STREAM_PARAMS = dict(
        mu_phi1cosphi2_prog=-0.38297458,
        mu_phi2_prog=-0.87059476,
        rv_prog=-109.48359169,
        dist_prog=21.8659734,
        phi2_prog=0.70106313,
        M_LMC=15,
    )

    # check if the necessary directories exist, if not, create them
    os.makedirs('orbits', exist_ok=True)
    os.makedirs('pre_impact', exist_ok=True)
    os.makedirs('final_stream', exist_ok=True)
    os.makedirs('final_coords', exist_ok=True)

    # Start simulation
    try:
        SH_x, SH_y, SH_z, SH_vx, SH_vy, SH_vz, _ = initial_stream.chi2_eval(
            tapproach=t_a, tmax=tmax, pid=int(pid), **DEFAULT_STREAM_PARAMS)
    except Exception as e:
        raise Exception(f"Error in initial_stream: {e}")
    try:
        x_sat, y_sat, z_sat, vx_sat, vy_sat, vz_sat = subhalo_orbit.chi2_eval(
            SH_x, SH_y, SH_z, SH_vx, SH_vy, SH_vz, r, phi, vphi, vz,
            tmax, t_a, phi_a, int(pid))
    except Exception as e:
        raise Exception(f"Error in subhalo_orbit: {e}")
    try:
        chi = stream_impact.chi2_eval(
            x_sat=x_sat, y_sat=y_sat, z_sat=z_sat, vx_sat=vx_sat, vy_sat=vy_sat,
            vz_sat=vz_sat, tmax=tmax, M_sat=M_sat, sr=rs_sat, pid=int(pid),
            **DEFAULT_STREAM_PARAMS)
    except Exception as e:
        raise Exception(f"Error in stream_impact: {e}")

    return chi

def process_stream(pid, params):
    """ Process the stream data into observables """

    os.makedirs('final_observables', exist_ok=True)

    # read in the stream data in RA/Dec coordinates
    data = np.genfromtxt(f'final_stream/final_stream_{pid}.txt')

    # transform to stream coordinates
    # TODO: this hack needs to be streamlined better
    R_PHI12_RADEC = np.array(
        [[0.83697865, 0.29481904, -0.4610298],
        [0.51616778, -0.70514011, 0.4861566],
        [0.18176238, 0.64487142, 0.74236331]]
    )
    phi1, phi2, dist, pm1, pm2, vr = obs_from_pos6d(
        data[:,:3], data[:,3:6], R_PHI12_RADEC)

    with h5py.File(f'final_observables/{pid}.hdf5', 'w') as f:
        f.create_dataset('phi1', data=phi1, compression='gzip')
        f.create_dataset('phi2', data=phi2, compression='gzip')
        f.create_dataset('dist', data=dist, compression='gzip')
        f.create_dataset('pm1', data=pm1, compression='gzip')
        f.create_dataset('pm2', data=pm2, compression='gzip')
        f.create_dataset('vr', data=vr, compression='gzip')
        f.attrs.update(params)

def simulate_dataset(config: config_dict.ConfigDict, workdir: str):
    """ Simulate a dataset of streams with subhalo impacts given the parameters """

    # 1. Preparing the working directory: don't change this
    # create working directories for the simulation
    logging.info(f'Creating working directory {workdir}')
    if os.path.exists(workdir) and (not config.overwrite):
        raise ValueError(f'Working directory {workdir} already exists')
    else:
        shutil.rmtree(workdir, ignore_errors=True)
    os.makedirs(workdir, exist_ok=True)

    # copy the binaries and potential files to the working directory
    logging.info('Copying binaries and potential files')
    shutil.copytree(config.binary_dir, os.path.join(workdir, 'binaries'))
    shutil.copytree(config.potential_dir, os.path.join(workdir, 'pot'))

    os.chdir(workdir)
    os.makedirs('orbits', exist_ok=True)
    os.makedirs('pre_impact', exist_ok=True)
    os.makedirs('final_stream', exist_ok=True)
    os.makedirs('final_coords', exist_ok=True)
    os.makedirs('final_observables', exist_ok=True)

    # 2. Reading in the Posterior file
    params_dict, num_params = read_posteriors(
        config.priors.filename,
        num_samples=config.priors.num_samples,
        quantiles=config.priors.quantiles,
        seed=config.seed
    )
    if config.priors.quantiles is None:
        batch_size = num_params // config.n_job
        start_idx = config.i_job * batch_size
        end_idx = start_idx + batch_size if config.i_job < config.n_job - 1 else num_params
    else:
        start_idx = 0
        end_idx = len(config.priors.quantiles)

    print(start_idx, end_idx)

    assert 1 == 2

    print(f'Run start-end: {start_idx}, {end_idx}')

    # 3. Go over the parameters one by one and simulate them
    all_params = {'chi2': []}  # all valid parameters
    for curr_pid in range(start_idx, end_idx):
        params = {k: params_dict[k][curr_pid] for k in params_dict}
        params['pid'] = int(curr_pid)
        try:
            # run simulation and process
            chi2 = run_single_sim(**params)
            process_stream(curr_pid, params)

            # add parameters and chi
            for key, val in params.items():
                if key not in all_params:
                    all_params[key] = []
                all_params[key].append(val)
            all_params['chi2'].append(chi2)
            print(all_params)

            # delete the pre-impact and final stream files
            if config.save_observables_only:
                os.remove('orbits/orbits_{:d}.txt'.format(curr_pid))
                os.remove('pre_impact/pre_impact_{:d}.txt'.format(curr_pid))
                os.remove('final_stream/final_stream_{:d}.txt'.format(curr_pid))
                os.remove('final_coords/final_coords_{:d}.txt'.format(curr_pid))

        except Exception as e:
            logging.warning(
                f'Failed to simulate stream {curr_pid} with parameters {params}')
            logging.warning(e)

    # save the simulation parameters as a CSV table
    logging.info('Saving simulation parameters...')
    all_params = pd.DataFrame(all_params)
    all_params.to_csv('labels.csv', index=False)

    # move the simulation results to the working directory
    logging.info('Finalizing simulation results...')
    if config.save_observables_only:
        logging.info('Saving only the observables. Deleting other files...')
        shutil.rmtree('orbits')
        shutil.rmtree('pre_impact')
        shutil.rmtree('final_stream')
        shutil.rmtree('final_coords')
        shutil.rmtree('binaries')
        shutil.rmtree('pot')

    return all_params


if __name__ == "__main__":
    FLAGS = flags.FLAGS
    config_flags.DEFINE_config_file(
        "config",
        None,
        "File path to the simulation configuration.",
        lock_config=True,
    )
    # Parse flags
    FLAGS(sys.argv)

    # Start training run
    t1 = time.time()
    simulate_dataset(config=FLAGS.config, workdir=FLAGS.config.workdir)
    t2 = time.time()

    logging.info(f'Simulation completed in {t2-t1:.2f} seconds')