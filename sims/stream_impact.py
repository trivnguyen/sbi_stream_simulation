
import os
import subprocess
from os.path import join

import astropy.units as u
import gal_uvw
import numpy
import numpy as np
import scipy
import scipy.optimize as sopt
from astropy.coordinates import SkyCoord
from astropy.io import fits
from rotation_matrix import phi12_rotmat, pmphi12

DEFAULT_BIN_PATH = "./binaries/stream_plummer.out"

# DEFINE CONSTANTS
M_NFW = 80.
q_NFW = 1.
rs_NFW = 16.
phi1_prog = 0.
phi2_prog = 0.
mu_phi2_prog = 0.
R_phi12_radec = np.array([[0.83697865, 0.29481904, -0.4610298],
                          [0.51616778, -0.70514011, 0.4861566],
                          [0.18176238, 0.64487142, 0.74236331]])
a_g = np.array([[-0.0548755604, +0.4941094279, -0.8676661490],
                [-0.8734370902, -0.4448296300, -0.1980763734],
                [-0.4838350155, 0.7469822445, +0.4559837762]])
G = 43007.105731706317
l_lmc = 280.4652*np.pi/180.
b_lmc = -32.8884*np.pi/180.

mcm = np.array([ 6.79343e+08,  2.82344e+00,  3.00000e-01,  0.00000e+00,
        0.00000e+00,  2.31801e+08,  2.95513e+00,  9.00000e-01,
        0.00000e+00,  0.00000e+00,  5.31319e+07,  7.00000e+00,
       -8.50000e-02,  4.00000e+00,  0.00000e+00,  2.17995e+09,
        1.50000e+00, -4.50000e-02,  1.20000e+01,  0.00000e+00,
        1.04699e+11,  5.00000e-01,  0.00000e+00,  1.80000e+00,
        7.50000e-02,  2.10000e+00,  1.57622e+07,  1.00000e+00,
        1.00000e+00,  3.00000e+00,  1.31445e+01,  0.00000e+00,
        8.40569e-03,  1.20149e-02,  7.28025e-03,  7.56561e-03,
       -3.34771e-03, -1.76854e-03,  1.92219e-03,  6.80997e+10,
        4.71197e+11,  6.38877e+11,  1.98064e+02,  8.95390e+11,
        1.96559e+01,  8.22844e+00,  2.33670e+02,  5.48871e+02,
        2.19013e-03,  1.60315e+03,  5.00000e+00])

#================================ATLAS DATA LOADED HERE, IGNORE FOR NOW========================================
#==============================================================================================================
fileout = open('pot/PJM17_{0}.Tpot'.format(0),'w')

print(4, file=fileout)
print('{0:.5e} {1:.5f} {2:.5f} {3} {4}'.format(mcm[0],mcm[1],mcm[2],mcm[3],mcm[4]), file=fileout)
print('{0:.5e} {1:.5f} {2:.5f} {3} {4}'.format(mcm[5],mcm[6],mcm[7],mcm[8],mcm[9]), file=fileout)
print('5.31319e+07 7 -0.085 4 0', file=fileout)
print('2.17995e+09 1.5 -0.045 12 0', file=fileout)

print(2, file=fileout)
print('{0:.5e} {1:.5f} {2:.5f} {3} {4} {5}'.format(mcm[20],mcm[21],mcm[22],mcm[23],mcm[24],mcm[25]), file=fileout)
print('{0:.5e} {1:.5f} {2:.5f} {3} {4} {5}'.format(mcm[26],mcm[27],mcm[28],mcm[29],mcm[30],mcm[31]), file=fileout)

Usun = mcm[32]*1000.
Vsun = mcm[33]*1000.
Wsun = mcm[34]*1000.
R0 = mcm[-6]
V0 = mcm[-5]

M200 = 4.*np.pi*mcm[26]*mcm[30]**3.*(np.log(1.+mcm[-9]/mcm[30])-mcm[-9]/(mcm[-9]+mcm[30]))/(1.e10)

c200 = mcm[-9]/mcm[30]
rs = mcm[30]

fileout.close()

def chi2_eval(
    mu_phi1cosphi2_prog, mu_phi2_prog, rv_prog,dist_prog, phi2_prog, M_LMC,
    x_sat, y_sat, z_sat, vx_sat, vy_sat, vz_sat, tmax, M_sat, sr, pid,
    bin_path=DEFAULT_BIN_PATH):

    # M_sat = 0.001
    print("Impact with scale radius: ",sr)

    mu_alpha_lmc = np.random.normal(1.91,0.*0.02)
    mu_delta_lmc = np.random.normal(0.229,0.*0.047)
    rv_lmc = np.random.normal(262.2,0.*3.4)
    dist_lmc = np.random.normal(49970.,0.*1126.)

    if M_LMC > 2.:
        rs_LMC = np.sqrt(G*M_LMC*8.7/91.7**2.)-8.7
    else:
        rs_LMC = np.sqrt(G*2.*8.7/91.7**2.)-8.7

    Mprog = 2.e-6

    lhood = chi2_worker(
        mu_phi1cosphi2_prog, mu_phi2_prog, rv_prog, dist_prog, phi2_prog,
        mu_alpha_lmc, mu_delta_lmc, rv_lmc, dist_lmc, M_LMC,rs_LMC, Mprog,
        tmax, x_sat,y_sat,z_sat,vx_sat,vy_sat,vz_sat,pid,M_sat, sr,
        bin_path=bin_path
    )
    return lhood

def chi2_worker(
    mu_phi1cosphi2_prog, mu_phi2_prog, rv_prog, dist_prog, phi2_prog,
    mu_alpha_lmc, mu_delta_lmc, rv_lmc, dist_lmc, M_LMC, rs_LMC, Mprog, tmax,
    x_sat, y_sat, z_sat, vx_sat, vy_sat, vz_sat, pid, M_sat, sr,
    bin_path=DEFAULT_BIN_PATH
):

    def loglikelihood(phi,phiM,errphiM):
        llh = np.sum(np.log((1/((2*np.pi*errphiM**2)**0.5)))-(((phiM-phi)**2)/(2*errphiM**2)))
        return llh

    def loglikelihoodQ(phi,phiM,errphiM):
        llh = np.sum(np.log((1/((2*np.pi*(errphiM**2+4))**0.5)))-(((phiM-phi)**2)/(2*(errphiM**2+4))))
        return llh

    vec_phi12_prog = np.array([
        np.cos(phi1_prog*np.pi/180.)*np.cos(phi2_prog*np.pi/180.),
        np.sin(phi1_prog*np.pi/180.)*np.cos(phi2_prog*np.pi/180.),
        np.sin(phi2_prog*np.pi/180.)
    ])

    vec_radec_prog = np.linalg.solve(R_phi12_radec,vec_phi12_prog)

    ra_prog = np.arctan2(vec_radec_prog[1],vec_radec_prog[0])*180./np.pi
    dec_prog = np.arcsin(vec_radec_prog[2]/np.linalg.norm(vec_radec_prog))*180./np.pi

    gc_prog = SkyCoord(ra=ra_prog*u.degree,dec=dec_prog*u.degree,frame='icrs')

    l_prog = np.array(gc_prog.galactic.l)
    b_prog = np.array(gc_prog.galactic.b)

    vlsr = np.array([Usun,Vsun+V0,Wsun])

    cos_phi1 = np.cos(phi1_prog * np.pi / 180.)
    sin_phi1 = np.sin(phi1_prog * np.pi / 180.)
    cos_phi2 = np.cos(phi2_prog * np.pi / 180.)
    sin_phi2 = np.sin(phi2_prog * np.pi / 180.)
    M_UVW_muphi12_prog = np.array([
        [cos_phi1 * cos_phi2, -sin_phi1, -cos_phi1 * sin_phi2],
        [sin_phi1 * cos_phi2, cos_phi1, -sin_phi1 * sin_phi2],
        [sin_phi2, 0. , cos_phi2]])

    k_mu = 4.74047

    uvw_stationary = -vlsr

    vec_vr_muphi1_muphi2_stationary = np.dot(
        M_UVW_muphi12_prog.T,np.dot(R_phi12_radec,np.dot(a_g,uvw_stationary)))
    vec_vr_muphi1_muphi2_stationary[0] = 0. # no correction for radial velocity, i want our radial velocity to be the los one

    vec_vr_muphi1_muphi2_prog = np.array([rv_prog,k_mu*dist_prog*mu_phi1cosphi2_prog,k_mu*dist_prog*mu_phi2_prog]) #+ vec_vr_muphi1_muphi2_stationary

    vx_prog,vy_prog,vz_prog = np.dot(a_g.T,np.dot(R_phi12_radec.T,np.dot(M_UVW_muphi12_prog,vec_vr_muphi1_muphi2_prog))) + vlsr

    x_prog, y_prog, z_prog = np.array([-R0,0.,0.])+dist_prog*np.array([np.cos(l_prog*np.pi/180.)*np.cos(b_prog*np.pi/180.),np.sin(l_prog*np.pi/180.)*np.cos(b_prog*np.pi/180.),np.sin(b_prog*np.pi/180.)])

    vx,vy,vz = vx_prog,vy_prog,vz_prog

    x,y,z = x_prog,y_prog,z_prog

    gc = SkyCoord(b=b_lmc*u.radian,l=l_lmc*u.radian,frame='galactic')

    x_lmc,y_lmc,z_lmc = np.array([-R0,0.,0.])+dist_lmc/1000.*np.array([np.cos(l_lmc)*np.cos(b_lmc),np.sin(l_lmc)*np.cos(b_lmc),np.sin(b_lmc)])

    vx_lmc,vy_lmc,vz_lmc = gal_uvw.gal_uvw(distance=dist_lmc,ra=np.array(gc.icrs.ra),dec=np.array(gc.icrs.dec),lsr=np.array([-Usun,0.,Wsun]),pmra=mu_alpha_lmc,pmdec=mu_delta_lmc,vrad=rv_lmc)

    vy_lmc += Vsun+V0
    vx_lmc = -vx_lmc

    # Run binary and read data
    cmd = '{28} {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} '\
        '{16} {17} {18} {19} {20} {21} {22} {23} {24} {25} {26}, {27}'.format(
            x, y, z, vx, vy, vz, M_LMC, rs_LMC, x_lmc, y_lmc, z_lmc,
            vx_lmc, vy_lmc, vz_lmc, Mprog, tmax, M200, rs, c200, x_sat, y_sat, z_sat,
            vx_sat, vy_sat, vz_sat, pid, M_sat, sr, bin_path
        )
    subprocess.check_call(cmd.split(' '))
    data = np.genfromtxt('final_stream/final_stream_{}.txt'.format(pid))

    pos = data[:,:3]
    vel = data[:,3:6]

    pos = pos - np.array([-R0,0.,0.])

    theta = np.arctan2(pos[:,1],pos[:,0])
    theta = np.mod(theta,2.*np.pi)
    theta[theta > np.pi] -= 2.*np.pi
    phi = np.arcsin(pos[:,2]/(pos[:,0]**2.+pos[:,1]**2.+pos[:,2]**2.)**0.5)

    l = 180./np.pi*theta
    b = 180./np.pi*phi

    gc = SkyCoord(l=l*u.degree,b=b*u.degree,frame='galactic')

    vec_radec = np.array([np.cos(gc.icrs.ra)*np.cos(gc.icrs.dec),np.sin(gc.icrs.ra)*np.cos(gc.icrs.dec),np.sin(gc.icrs.dec)])

    vec_phi12 = np.dot(R_phi12_radec,vec_radec).T

    phi1 = np.arctan2(vec_phi12[:,1],vec_phi12[:,0])*180./np.pi
    phi2 = np.arcsin(vec_phi12[:,2])*180./np.pi

    vel -= vlsr

    #ATLAS
    # alpha = np.array(gc.icrs.ra)*np.pi/180.
    # delta = np.array(gc.icrs.dec)*np.pi/180.

    r_stream = sum(pos[:,axis]**2. for axis in range(3))**0.5

    R_phi12_a_g = np.dot(R_phi12_radec,a_g)

    vr_stream = sum((np.cos(phi1*np.pi/180.)*np.cos(phi2*np.pi/180.)*R_phi12_a_g[0,axis]+np.sin(phi1*np.pi/180.)*np.cos(phi2*np.pi/180.)*R_phi12_a_g[1,axis]+np.sin(phi2*np.pi/180.)*R_phi12_a_g[2,axis])*vel[:,axis] for axis in range(3))

    mu_phi1_cos_phi2_stream = 1./(k_mu*r_stream)*sum( (-np.sin(phi1*np.pi/180.)*R_phi12_a_g[0,axis]+np.cos(phi1*np.pi/180.)*R_phi12_a_g[1,axis])*vel[:,axis] for axis in range(3))

    mu_phi2_stream = 1./(k_mu*r_stream)*sum( (-np.cos(phi1*np.pi/180.)*np.sin(phi2*np.pi/180.)*R_phi12_a_g[0,axis] - np.sin(phi1*np.pi/180.)*np.sin(phi2*np.pi/180.)*R_phi12_a_g[1,axis] + np.cos(phi2*np.pi/180.)*R_phi12_a_g[2,axis])*vel[:,axis] for axis in range(3))


    return 0

'''chi2_eval(-0.38297458,   -0.87059476, -109.48359169,   21.8659734 , 0.70106313,15,-20.010377,   -1.965021,  -10.462969, -107.437575, -189.452041,  -49.937659,2)'''


'''chi2_eval(-0.38297458,   -0.87059476, -109.48359169,   21.8659734 , 0.70106313,15,-8.561333,  14.537095, -17.848972,  97.351609, 189.817463, 142.413466,101)

chi2_eval(-0.38297458,   -0.87059476, -109.48359169,   21.8659734 , 0.70106313,15,-15.795487, -32.932709, -29.688411,  41.280606,  92.595954,   2.978653,102)

chi2_eval(-0.38297458,   -0.87059476, -109.48359169,   21.8659734 , 0.70106313,15,  5.388922,  -4.531722, -23.327607,  58.941398 ,193.300301,  93.196695,103)

chi2_eval(-0.38297458,   -0.87059476, -109.48359169,   21.8659734 , 0.70106313,15,-35.807043, -10.809626, -27.906099, -17.374411, 144.187285,  19.687321,104)

chi2_eval(-0.38297458,   -0.87059476, -109.48359169,   21.8659734 , 0.70106313,15,-6.295656,   1.005379,  -0.450091, 252.385492, 286.896465, 189.446385,105)

chi2_eval(-0.38297458,   -0.87059476, -109.48359169,   21.8659734 , 0.70106313,15,-14.122377, -11.27479,  -50.172046,  40.3789,   133.469086, -39.612097,106)'''
if __name__=='__main__':

   #chi2_eval(mu_phi1cosphi2_prog, mu_phi2_prog, rv_prog,dist_prog,phi2_prog,M_LMC,x_sat,y_sat,z_sat,vx_sat,vy_sat,vz_sat,tmax,M_sat,pid):
    chi2_eval(-0.38297458,   -0.87059476, -109.48359169,   21.8659734 , 0.70106313,15,-27.150089,-14.180522,-65.251477,-19.181674,115.954064,-109.696156,3,0.001,0)

    # chi2_eval(-0.38297458,   -0.87059476, -109.48359169,   21.8659734 , 0.70106313,15,-13.855882,-11.527033,-42.728755,48.538603,136.61262,-9.013355,0.00075,100555)
    # chi2_eval(-0.38297458,   -0.87059476, -109.48359169,   21.8659734 , 0.70106313,15,-4.703102,-20.535064,-23.575307,69.868852,147.045624,39.728405,0.0005,10066)

    # chi2_eval(-0.38297458,   -0.87059476, -109.48359169,   21.8659734 , 0.70106313,15,-8.089494,-16.237244,-31.54951,70.734234,138.122854,27.597337,0.002,10055)

    # chi2_eval(-0.38297458,   -0.87059476, -109.48359169,   21.8659734 , 0.70106313,15,-15.134876,-9.588946,-45.217458,42.73426,139.052821,-16.777925,0.001,1003)
    # chi2_eval(-0.38297458,   -0.87059476, -109.48359169,   21.8659734 , 0.70106313,15,-15.134876,-9.588946,-45.217458,42.73426,139.052821,-16.777925,0.00075,10044)

    # chi2_eval(-0.38297458,   -0.87059476, -109.48359169,   21.8659734 , 0.70106313,15,-8.33438,15.062839,-17.556942,98.015064,189.15073,146.646894,0.001,101)

    # chi2_eval(-0.38297458,   -0.87059476, -109.48359169,   21.8659734 , 0.70106313,15,-15.402457,-32.018567,-29.379787,44.254697,95.956652,7.813829,0.001,102)

    # chi2_eval(-0.38297458,   -0.87059476, -109.48359169,   21.8659734 , 0.70106313,15,5.649581,-3.811265,-22.887055,59.186669,195.009433,100.899814,0.001,103)

    # chi2_eval(-0.38297458,   -0.87059476, -109.48359169,   21.8659734 , 0.70106313,15,-35.449946,-9.994334,-27.586293,-14.539456,146.940119,24.801512,0.001,104)

    # chi2_eval(-0.38297458,   -0.87059476, -109.48359169,   21.8659734 , 0.70106313,15,-5.444992,1.889605,0.239211,278.448404,280.564033,195.004059,0.001,105)

    # chi2_eval(-0.38297458,   -0.87059476, -109.48359169,   21.8659734 , 0.70106313,15,-13.781209,-10.508278,-49.918173,42.332814,135.641759,-35.040883,0.001,106)

    # print('module 3: this should not appear')
