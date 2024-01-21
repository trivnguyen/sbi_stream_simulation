# from subhalo_impact import chi_eval
from multiprocessing import Process
from multiprocessing import Pool
import itertools

import initial_stream
import subhalo_orbit
import stream_impact
from rotation_matrix import obs_from_pos6d
import numpy as np
import h5py
import os

'''
0-1699 M, vz
1700-2199 M, vz, r
2200-8499 M (0.1), vz
8500-14499 M(0.1), vz (10<|vz|<200), r (0.1-1)

14500-14520 test Mlog=-5
14520-14540 test Mlog=-1
14540-14560 test |vz| = 5
14560-14580 test |vz| = 300
14580-14600 test r = 0.01
14600-14620 test r = 1
14620-14640 test t_a = 0
14640-14660 test t_a = 4
14660-14670 test phi = -180
14670-14680 test phi = 0
14680-14690 test phi = 180
14690-14700 test phi_a = -30
14700-14710 test phi_a = 0
14710-14720 test phi_a = 30
14720-14740 test rs = 0.01
14740-14760 test rs = 5
14760-14780 test vphi = 5 ?
14780-14800 test vphi = 300
'''

def run_sims(nsims=1000):
    pool = Pool()

    '''
    logM, vz, r, ta, phi, phia, rs, vphi
    '''
    logM_array = np.random.uniform(-5,-1,nsims)
    vz_factor = np.random.choice([-1, 1], size = nsims)
    vz_array = np.random.uniform(5,300,nsims)
    vz_array = vz_factor * vz_array
    r_array = np.random.uniform(0.01, 1, nsims)
    t_a_array = np.random.uniform(0,4,nsims)
    phi_array = np.random.uniform(-180, 180, nsims)
    phi_a_array = np.random.uniform(-25, 25, nsims)
    rs_array = np.random.uniform(0.01, 5, nsims)
    v_phi_factor = np.random.choice([-1, 1], size = nsims)
    v_phi_array = np.random.uniform(5, 6, nsims)
    v_phi_array = v_phi_array * v_phi_factor
    n = 14763
    pid_array = list(range(n, n + nsims))

    
    '''
    logM_array = np.random.uniform(-5,-1,nsims)
    vz_factor = np.random.choice([-1, 1], size = nsims)
    vz_array = np.random.uniform(5,300,nsims)
    vz_array = vz_factor * vz_array
    r_array = np.random.uniform(0.01, 1, nsims)
    t_a_array = np.random.uniform(0,4,nsims)
    phi_array = np.random.uniform(-180, 180, nsims)
    phi_a_array = np.random.uniform(-30, 30, nsims)
    rs_array = np.random.uniform(0.01, 5, nsims)
    v_phi_factor = np.random.choice([-1, 1], size = nsims)
    v_phi_array = np.random.uniform(5, 300, nsims)
    v_phi_array = v_phi_array * v_phi_factor
    '''

    pool.map(run_sim, zip(logM_array, vz_array, r_array, t_a_array, phi_array, phi_a_array, rs_array, v_phi_array, pid_array))
    
def run_sim(params):
    logM_sat, vz, r, t_a, phi, phi_a, rs_sat, vphi, pid = params
    M_sat = 10**logM_sat
    # rs_sat = 1.05 * (M_sat*10*10)**0.5 
    
    # r = 0.2 # impact parameter in kpc (distance from stream to sat)
    # phi = 250 # angle around stream in dec
    # vphi = 35  # velocity around stream in km/s
    # t_a = 0.2 # time since interaction in Gyr
    # phi_a = -4 # interaction point along stream in deg, phi=0 is progenitor location (try -20 to 10)
    tmax=4 # how long stream disrupts in Gyr

    print()
    print()
    print('Running logM = %.3f, vz = %.3f' %(M_sat, vz))
    simulate_stream(r,phi,vphi,vz,M_sat,tmax,t_a,phi_a,rs_sat,pid)


def simulate_stream(r, phi, vphi, vz, M_sat, tmax, t_a, phi_a, rs_sat, pid):
    SH_x, SH_y, SH_z, SH_vx, SH_vy, SH_vz, dunno = initial_stream.chi2_eval(-0.38297458,   -0.87059476, -109.48359169,   21.8659734 , 0.70106313,15,t_a,tmax,int(pid))
    sat_x, sat_y, sat_z, sat_vx, sat_vy, sat_vz = subhalo_orbit.chi2_eval(SH_x, SH_y, SH_z, SH_vx, SH_vy, SH_vz,r,phi,vphi,vz,tmax,t_a,phi_a,int(pid))
    chi = stream_impact.chi2_eval(-0.38297458,   -0.87059476, -109.48359169,   21.8659734 , 0.70106313,15, sat_x, sat_y, sat_z, sat_vx, sat_vy, sat_vz, tmax,M_sat,rs_sat,int(pid))

    # save observables as a new file after simulating the stream

    save_observables_hdf5(r, phi, vphi, vz, M_sat, tmax, t_a, phi_a, rs_sat, pid)

    # remove waste files
    os.remove('orbits/orbit_%i.txt' %pid)
    os.remove('pre_impact/pre_impact_%i.txt' %pid)
    os.remove('final_coords/final_coords_%i.txt' %pid)
    os.remove('final_stream/final_stream_%i.txt' %pid)


R_phi12_radec = np.array([[0.83697865, 0.29481904, -0.4610298], 
                          [0.51616778, -0.70514011, 0.4861566], 
                          [0.18176238, 0.64487142, 0.74236331]])

def save_observables_hdf5(r, phi, vphi, vz, M_sat, tmax, t_a, phi_a, rs_sat, pid):
    data = np.genfromtxt(f'final_stream/final_stream_{pid}.txt')
    phi1,phi2,dist,pm1,pm2,vr = obs_from_pos6d(data[:,:3],data[:,3:6],R_phi12_radec)
    phi1 = phi1.astype('float32')
    phi2 = phi2.astype('float32')
    dist = dist.astype('float32')
    pm1 = pm1.astype('float32')
    pm2 = pm2.astype('float32')
    vr = vr.astype('float32')

    parameters = np.array([r, phi, vphi, vz, M_sat, tmax, t_a, phi_a, rs_sat])
    parameters = parameters.astype('float32')

    f = h5py.File(f'observables/sim{pid}.hdf5', 'w')
    f.create_dataset('parameters', data = parameters, compression = 'gzip')
    f.create_dataset('phi1', data = phi1, compression = 'gzip')
    f.create_dataset('phi2', data = phi2, compression = 'gzip')
    f.create_dataset('dist', data = dist, compression = 'gzip')
    f.create_dataset('pm1', data = pm1, compression = 'gzip')
    f.create_dataset('pm2', data = pm2, compression = 'gzip')
    f.create_dataset('vr', data = vr, compression = 'gzip')
    f.close()

def read_observables_hdf5(hf):
    f = h5py.File(f'observables/{hf}', 'r')
    parameters = np.array(f.get('parameters'))
    phi1 = np.array(f.get('phi1'))
    phi2 = np.array(f.get('phi2'))
    dist = np.array(f.get('dist'))
    pm1 = np.array(f.get('pm1'))
    pm2 = np.array(f.get('pm2'))
    vr = np.array(f.get('vr'))
    f.close()
    return parameters, phi1, phi2, dist, pm1, pm2, vr

def save_pid_txt(pid):
    filename = f'observables/sim{pid}.hdf5'
    with h5py.File(f'{filename}', "r") as hdf5_file:
            data = np.array(hdf5_file.get('parameters'))
            r = data[0]
            phi = data[1]
            vphi = data[2]
            vz = data[3]
            M_sat = data[4]
            tmax = data[5]
            t_a = data[6]
            phi_a = data[7]
            rs_sat = data[8]
            hdf5_file.close()

            with open("pid.txt", "a") as f:
                f.write(f"{pid}, {r}, {phi}, {vphi}, {vz}, {M_sat}, {tmax}, {t_a}, {phi_a}, {rs_sat}\n")
                f.close()

# '''re-sorttt_everying_in_pidtxt_ahhhh''':
# import temp
# import os
# import glob
# os.chdir(r'/home/rutong/stellar_stream/SBI-streameeeee/sims/observables')
# my_files = glob.glob('*.hdf5')
# os.chdir(r'/home/rutong/stellar_stream/SBI-streameeeee/sims')
# pid_list = [int(f[3:-5]) for f in my_files]
# list.sort(pid_list)
# for pid in pid_list:
#     temp.save_pid_txt(pid)


'''
def load_parameters_txt(??):
    data = np.loadtxt("pid.txt", delimiter=",", skiprows=1)
    second_column = data[:, 1]
'''

if __name__ == '__main__':
    run_sims(nsims=1)



