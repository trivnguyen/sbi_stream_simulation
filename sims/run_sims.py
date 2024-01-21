import itertools

import initial_stream
import subhalo_orbit
import stream_impact
from rotation_matrix import obs_from_pos6d
import numpy as np
import os
import h5py


# def run_sims(nsims=1000):
#     pool = Pool()
#     logM_array = np.random.uniform(-5,1,nsims)
#     vz_array = np.random.uniform(-50,0,nsims)
#     pool.map(run_sim, zip(logM_array, vz_array))

def run_sim(params, pid):
    logM_sat, vz = params
    M_sat = 10**logM_sat
    # M_sat = 10**float(sys.argv[1])
    # vz = 10**float(sys.argv[2])
    rs_sat = 1.05 * (M_sat*10*10)**0.5

    print(pid)

    # pid = calculate_pid(logM_sat, vz)

    r = 0.2 # impact parameter in kpc (distance from stream to sat)
    phi = 250 # angle around stream in dec
    vphi = 35  # velocity around stream in km/s
    # vz = -10  # velocity along stream in km/s
    # M_sat = 0.001 # mass of subhalo in 1e10 Msun
    t_a = 0.2 # time since interaction in Gyr
    phi_a = -4 # interaction point along stream in deg, phi=0 is progenitor location (try -20 to 10)
    # rs_sat = 0.3 # scale radius of subhalo in kpc, can be adjusted along with M_sat using equation 15 in erkal et al. 2015
    # pid=0 # index included in saved filenames
    tmax=4 # how long stream disrupts in Gyr

    print()
    print()
    print('Running logM = %.3f, vz = %.3f' %(M_sat, vz))
    print()
    print()
    simulate_stream(r,phi,vphi,vz,M_sat,tmax,t_a,phi_a,rs_sat,pid)


def simulate_stream(r, phi, vphi, vz, M_sat, tmax, t_a, phi_a, rs_sat, pid):
    print('0')
    SH_x, SH_y, SH_z, SH_vx, SH_vy, SH_vz, dunno = initial_stream.chi2_eval(-0.38297458,   -0.87059476, -109.48359169,   21.8659734 , 0.70106313,15,t_a,tmax,int(pid))
    print('1')
    sat_x, sat_y, sat_z, sat_vx, sat_vy, sat_vz = subhalo_orbit.chi2_eval(SH_x, SH_y, SH_z, SH_vx, SH_vy, SH_vz,r,phi,vphi,vz,tmax,t_a,phi_a,int(pid))
    print('2')
    chi = stream_impact.chi2_eval(-0.38297458,   -0.87059476, -109.48359169,   21.8659734 , 0.70106313,15, sat_x, sat_y, sat_z, sat_vx, sat_vy, sat_vz, tmax,M_sat,rs_sat,int(pid))

    os.remove('orbits/orbit_%i.txt' %pid)
    os.remove('pre_impact/pre_impact_%i.txt' %pid)
    os.remove('final_coords/final_coords_%i.txt' %pid)

    # save observables as a new file after simulating the stream

    # save_observables_hdf5(pid)
    save_observables_hdf5(pid)


R_phi12_radec = np.array([[0.83697865, 0.29481904, -0.4610298],
                          [0.51616778, -0.70514011, 0.4861566],
                          [0.18176238, 0.64487142, 0.74236331]])

def save_observables_txt(pid):
    data = np.genfromtxt(f'final_stream/final_stream_{pid}.txt')
    phi1,phi2,dist,pm1,pm2,vr = obs_from_pos6d(data[:,:3],data[:,3:6],R_phi12_radec)
    observables = np.vstack((phi1, phi2, dist, pm1, pm2, vr))
    np.savetxt(f'final_observables/{pid}.txt', observables)

def save_observables_hdf5(pid):
    data = np.genfromtxt(f'final_stream/final_stream_{pid}.txt')
    phi1,phi2,dist,pm1,pm2,vr = obs_from_pos6d(data[:,:3],data[:,3:6],R_phi12_radec)
    phi1 = phi1.astype('float32')
    phi2 = phi2.astype('float32')
    dist = dist.astype('float32')
    pm1 = pm1.astype('float32')
    pm2 = pm2.astype('float32')
    vr = vr.astype('float32')

    f = h5py.File(f'final_observables/{pid}.hdf5', 'w')
    f.create_dataset('phi1', data = phi1, compression = 'gzip')
    f.create_dataset('phi2', data = phi2, compression = 'gzip')
    f.create_dataset('dist', data = dist, compression = 'gzip')
    f.create_dataset('pm1', data = pm1, compression = 'gzip')
    f.create_dataset('pm2', data = pm2, compression = 'gzip')
    f.create_dataset('vr', data = vr, compression = 'gzip')
    f.close()

def read_observables_txt(txtfile):
    data = np.genfromtxt(txtfile)
    return data[0], data[1], data[2], data[3], data[4], data[5]

def read_observables_hdf5(hdf5file):
    f = h5py.File(hdf5file, 'r')
    phi1 = np.array(f.get('phi1'))
    phi2 = np.array(f.get('phi1'))
    dist = np.array(f.get('dist'))
    pm1 = np.array(f.get('phi1'))
    pm2 = np.array(f.get('phi1'))
    vr = np.array(f.get('vr'))
    f.close()
    return phi1, phi2, dist, pm1, pm2, vr

if __name__ == '__main__':
    nsims = 1000
    logM_array = np.random.uniform(-5, 1, nsims)
    vz_array = np.random.uniform(-50, 0, nsims)

    for i in range(0, 5):
        try:
            pid = i
            run_sim([logM_array[i], vz_array[i]], pid)
            print(i)
        except Exception as e:
            print(e)