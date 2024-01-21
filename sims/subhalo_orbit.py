
import os
import subprocess
from os.path import join

import astropy.units as u
import astropy.io.fits as afits
import csv
import numpy as np
import scipy.optimize as sopt
import scipy.interpolate as sinterp
import scipy
import gal_uvw
from scipy.interpolate import UnivariateSpline
from rotation_matrix import phi12_rotmat, pmphi12
from astropy.io import fits
from astropy.coordinates import SkyCoord

DEFAULT_BIN_PATH = './binaries/orbit.out'

# DEFAULT PARAMETERS
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

def impact_params(pid,tmax,tapproach,phiapproach,r_o,phi_o,vphi_o,vz_oi):

    if np.around(r_o,decimals=5)==0:
        direc = 1
    else:
        direc = r_o/abs(r_o)

    # print('r_o, direc:',r_o,direc)
    x_o = abs(r_o)*np.cos(np.radians(phi_o))
    y_o = 0
    z_o = abs(r_o)*np.sin(np.radians(phi_o))
    vx_o = -vphi_o*np.sin(np.radians(phi_o))*direc
    vy_o = vz_oi
    vz_o = vphi_o*np.cos(np.radians(phi_o))*direc

    # pid = 0
    # tmax = 3
    # tapproach = 0.25
    phiapproach = np.radians(-phiapproach)

    orbit_dat = np.genfromtxt(join('orbits','orbit_{0}.txt'.format(pid)))

    stream_pos = orbit_dat[:,2:5]
    stream_vel = orbit_dat[:,5:8]
    mw_pos = orbit_dat[:,14:17]
    mw_vel = orbit_dat[:,17:20]
    t = orbit_dat[:,1]/1.0227
    t = t.round(4)

    prog_pos = stream_pos-mw_pos
    prog_vel = -(stream_vel-mw_vel)

    #Now get positions and velocities at t_approach


    #TODO: make approach 2 indices either side, with a third value for the time difference?
    #np.array([np.interp(0.5,[0,1e-4],[vect1[0],vect2[0]])
            # ,np.interp(0.5,[0,1e-4],[vect1[1],vect2[1]])
            # ,np.interp(0.5,[0,1e-4],[vect1[2],vect2[2]])])

    #np.array([np.interp(t_diff,[0,1e-4],[prog_pos[target1,0],prog_pos[target2,0]])
            # ,np.interp(t_diff,[0,1e-4],[prog_pos[target1,1],prog_pos[target2,1]])
            # ,np.interp(t_diff,[0,1e-4],[prog_pos[target1,2],prog_pos[target2,2]])])

    targetapproach = tmax-tapproach
    target1 = t.size - np.searchsorted(t[::-1], targetapproach, side = "right") #this chooses the lower end t
    target2 = target1-1
    t_diff = targetapproach-t[target1]

    prog_pos = np.array([np.interp(t_diff,[0,1e-4],[prog_pos[target1,0],prog_pos[target2,0]])
                        ,np.interp(t_diff,[0,1e-4],[prog_pos[target1,1],prog_pos[target2,1]])
                        ,np.interp(t_diff,[0,1e-4],[prog_pos[target1,2],prog_pos[target2,2]])])

    prog_vel = np.array([np.interp(t_diff,[0,1e-4],[prog_vel[target1,0],prog_vel[target2,0]])
                        ,np.interp(t_diff,[0,1e-4],[prog_vel[target1,1],prog_vel[target2,1]])
                        ,np.interp(t_diff,[0,1e-4],[prog_vel[target1,2],prog_vel[target2,2]])])

    mw_pos_approach = np.array([np.interp(t_diff,[0,1e-4],[mw_pos[target1,0],mw_pos[target2,0]])
                               ,np.interp(t_diff,[0,1e-4],[mw_pos[target1,1],mw_pos[target2,1]])
                               ,np.interp(t_diff,[0,1e-4],[mw_pos[target1,2],mw_pos[target2,2]])])

    mw_vel_approach = np.array([np.interp(t_diff,[0,1e-4],[mw_vel[target1,0],mw_vel[target2,0]])
                               ,np.interp(t_diff,[0,1e-4],[mw_vel[target1,1],mw_vel[target2,1]])
                               ,np.interp(t_diff,[0,1e-4],[mw_vel[target1,2],mw_vel[target2,2]])])

    prog_L = np.cross(prog_pos,prog_vel)

    #First, rotate to x-z plane
    rtheta1 = ((prog_L[1]*prog_L[2])/abs((prog_L[1]*prog_L[2])))*np.arctan2(prog_L[1],prog_L[2])

    rm1 = np.array(([1,      0        ,       0        ],
                    [0,np.cos(rtheta1),-np.sin(rtheta1)],
                    [0,np.sin(rtheta1), np.cos(rtheta1)]))

    prog_L1 = np.matmul(rm1,prog_L)

    #Now rotate to z
    rtheta2 = ((prog_L[0]*prog_L[2])/abs((prog_L[0]*prog_L[2])))*np.arctan2(prog_L1[0],prog_L1[2])

    rm2 = np.array(([ np.cos(rtheta2),0,np.sin(rtheta2)],
                    [       0        ,1,      0        ],
                    [-np.sin(rtheta2),0,np.cos(rtheta2)]))

    #This rotation matrix can be used on the stream to move it mostly into the x-y plane
    rm3 = np.matmul(rm2,rm1)

    data_approach = np.genfromtxt(join('pre_impact','pre_impact_{0}.txt'.format(pid)))

    stream_pos = data_approach[:,0:3]-mw_pos_approach
    stream_vel = data_approach[:,3:6]+mw_vel_approach

    prog_pos2 = rm3.dot(prog_pos)
    stream_pos2 = rm3.dot(stream_pos.T).T
    stream_vel2 = rm3.dot(stream_vel.T).T

    phi_prog = -np.arctan2(prog_pos2[1],prog_pos2[0])

    # rmG = np.array(([np.cos(phi_prog),-np.sin(phi_prog),0],
    #                 [np.sin(phi_prog), np.cos(phi_prog),0],
    #                 [              0             ,                0            ,1]))
    rmG = np.array(([np.cos(phi_prog+phiapproach),-np.sin(phi_prog+phiapproach),0],
                    [np.sin(phi_prog+phiapproach), np.cos(phi_prog+phiapproach),0],
                    [              0             ,                0            ,1]))
    stream_posG = rmG.dot(stream_pos2.T).T
    stream_velG = rmG.dot(stream_vel2.T).T
    # print("Stream Velocities: ",stream_velG)
    stream_phi = np.arctan2(stream_posG[:,1],stream_posG[:,0])
    # print(stream_phi)

    vh_1 = np.radians(-0.25)
    vh_2 = np.radians(0.25)
    vh_stars = np.where(((stream_phi>=vh_1)&(stream_phi<=vh_2)))[0]
    # print("Stars in slice: ",vh_stars.shape)
    stream_posV = stream_posG[vh_stars]
    stream_velV = stream_velG[vh_stars]
    vh_x = np.mean(stream_posV[:,0])
    vh_y = np.mean(stream_posV[:,1])
    vh_z = np.mean(stream_posV[:,2])
    vh_vx = np.mean(stream_velV[:,0])
    vh_vy = np.mean(stream_velV[:,1])
    vh_vz = np.mean(stream_velV[:,2])

    #The impact parameters in the rotated frame
    impact_posG = np.array([vh_x+x_o,vh_y+y_o,vh_z+z_o])
    impact_velG = np.array([vh_vx+vx_o,vh_vy+vy_o,vh_vz+vz_o])



    impact_pos = np.linalg.inv(rmG).dot(impact_posG)
    impact_vel = np.linalg.inv(rmG).dot(impact_velG)

    impact_pos_o = np.linalg.inv(rm3).dot(impact_pos) + mw_pos_approach
    impact_vel_o = np.linalg.inv(rm3).dot(impact_vel) - mw_vel_approach

    #temporarily return b_out here
    # print('returning the following positions: ',impact_posG)
    # print('returning the following velocities: ',impact_velG)
    # print('Function using the following params: ',[pid,tmax,tapproach,phiapproach,r_o,phi_o,vphi_o,vz_oi])

    # print('vx',-vphi_o*np.sin(np.radians(phi_o))*r_o/abs(r_o))
    # print('vz',vphi_o*np.cos(np.radians(phi_o))*r_o/abs(r_o))

    # print('TESTING HERE========================', impact_pos_o[0], impact_pos_o[1], impact_pos_o[2], impact_vel_o[0], impact_vel_o[1], impact_vel_o[2])

    return impact_pos_o[0], impact_pos_o[1], impact_pos_o[2], impact_vel_o[0], impact_vel_o[1], impact_vel_o[2]
    # return brange, fullset


def chi2_eval(
    x_s, y_s, z_s, vx_s, vy_s, vz_s, r_o, phi_o, vphi_o, vz_o,
    tmax, tapproach, phiapproach, pid, bin_path=DEFAULT_BIN_PATH
):

    print('====================================================================================')
    print('Running orbit with r={},phi={},vphi={},vz={},tapproach={},phiapproach={},pid={}'.format(r_o,phi_o,vphi_o,vz_o,tapproach,phiapproach,pid))
    x_imp, y_imp, z_imp, vx_imp, vy_imp, vz_imp = impact_params(
        pid, tmax, tapproach, phiapproach, r_o,phi_o,vphi_o,vz_o)


    M_LMC = 15.

    mu_alpha_lmc = np.random.normal(1.91,0.*0.02)
    mu_delta_lmc = np.random.normal(0.229,0.*0.047)
    rv_lmc = np.random.normal(262.2,0.*3.4)
    dist_lmc = np.random.normal(49970.,0.*1126.)

    if M_LMC > 2.:
        rs_LMC = np.sqrt(G*M_LMC*8.7/91.7**2.)-8.7
    else:
        rs_LMC = np.sqrt(G*2.*8.7/91.7**2.)-8.7
    Mprog = 0.

    lhood = chi2_worker(
        x_s, y_s, z_s, vx_s, vy_s, vz_s, x_imp, y_imp, z_imp,
        vx_imp, vy_imp, vz_imp, mu_alpha_lmc, mu_delta_lmc, rv_lmc, dist_lmc,
        M_LMC, rs_LMC, Mprog, tmax,tapproach, pid,
        bin_path=bin_path
    )
    return lhood

def chi2_worker(
    x_s, y_s, z_s, vx_s, vy_s, vz_s, x_o, y_o, z_o, vx_o, vy_o, vz_o,
    mu_alpha_lmc, mu_delta_lmc, rv_lmc, dist_lmc, M_LMC, rs_LMC,
    Mprog, tmax, tapproach, pid, bin_path=DEFAULT_BIN_PATH
):

    #print('chi2 worker started')

    vlsr = np.array([Usun,Vsun+V0,Wsun])

    k_mu = 4.74047

    x=x_s
    y=y_s
    z=z_s
    vx = vx_s
    vy = vy_s
    vz = vz_s

    x_sat = x_o
    y_sat = y_o
    z_sat = z_o
    vx_sat = vx_o
    vy_sat = vy_o
    vz_sat = vz_o

    gc = SkyCoord(b=b_lmc*u.radian,l=l_lmc*u.radian,frame='galactic')

    x_lmc,y_lmc,z_lmc = np.array([-R0,0.,0.])+dist_lmc/1000.*np.array([np.cos(l_lmc)*np.cos(b_lmc),np.sin(l_lmc)*np.cos(b_lmc),np.sin(b_lmc)])

    vx_lmc,vy_lmc,vz_lmc = gal_uvw.gal_uvw(distance=dist_lmc,ra=np.array(gc.icrs.ra),dec=np.array(gc.icrs.dec),lsr=np.array([-Usun,0.,Wsun]),pmra=mu_alpha_lmc,pmdec=mu_delta_lmc,vrad=rv_lmc)

    vy_lmc += Vsun+V0
    vx_lmc = -vx_lmc

    # Run binaries
    cmd = '{26} {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} '\
        '{13} {14} {15} {16} {17} {18} {19} {20} {21} {22} {23} {24} {25}'.format(
            5, 5, 5, 1, 1, 1, M_LMC, rs_LMC, x_lmc, y_lmc, z_lmc, vx_lmc, vy_lmc, vz_lmc,
            Mprog, tapproach, M200, rs, c200, x_sat, y_sat, z_sat, vx_sat, vy_sat, vz_sat, pid,
            bin_path
        )
    subprocess.check_call(cmd.split(' '))
    data = np.genfromtxt('final_coords/final_coords_{}.txt'.format(pid))


    return data[0],data[1],data[2],data[3],data[4],data[5]


if __name__=='__main__':
    # plt.clf()
    sub_ics = np.array([-12.039907,-6.359913,-24.642776,71.477487,176.178943,73.943191,0.493398])

                                                                            #  r,phi,vphi,vz,tmax,tapproach,phiapproach,pid
    chi2_eval(sub_ics[0],sub_ics[1],sub_ics[2],sub_ics[3],sub_ics[4],sub_ics[5],0,270,50,100,4,6.46822173e-01,-3.5,10000)

    # chi2_eval(sub_ics[0],sub_ics[1],sub_ics[2],sub_ics[3],sub_ics[4],sub_ics[5],0,0,0,0,-100,0,0.25,2)

    # chi2_eval(sub_ics[0],sub_ics[1],sub_ics[2],sub_ics[3],sub_ics[4],sub_ics[5],0,0,0,100,0,0,0.25,3)

    # chi2_eval(sub_ics[0],sub_ics[1],sub_ics[2],sub_ics[3],sub_ics[4],sub_ics[5],0,0,0,-100,0,0,0.25,4)

    # chi2_eval(sub_ics[0],sub_ics[1],sub_ics[2],sub_ics[3],sub_ics[4],sub_ics[5],0,0,0,0,0,100,0.25,5)

    # chi2_eval(sub_ics[0],sub_ics[1],sub_ics[2],sub_ics[3],sub_ics[4],sub_ics[5],0,0,0,0,0,-100,0.25,6)
