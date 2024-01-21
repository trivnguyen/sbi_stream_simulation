import numpy as np


def get_potential(idx=9):
    mcmillan = np.genfromtxt('ProductionRunBig2_Nora_10.tab')
    mcm = mcmillan[idx]

    fileout = open('PJM17_{0}.Tpot'.format(0), 'w')

    print(4, file=fileout)
    print('{0:.5e} {1:.5f} {2:.5f} {3} {4}'.format(mcm[0], mcm[1], mcm[2], mcm[3], mcm[4]), file=fileout)
    print('{0:.5e} {1:.5f} {2:.5f} {3} {4}'.format(mcm[5], mcm[6], mcm[7], mcm[8], mcm[9]), file=fileout)
    print('5.31319e+07 7 -0.085 4 0', file=fileout)
    print('2.17995e+09 1.5 -0.045 12 0', file=fileout)

    print(2, file=fileout)
    print('{0:.5e} {1:.5f} {2:.5f} {3} {4} {5}'.format(mcm[20], mcm[21], mcm[22], mcm[23], mcm[24], mcm[25]), file=fileout)
    print('{0:.5e} {1:.5f} {2:.5f} {3} {4} {5}'.format(mcm[26], mcm[27], mcm[28], mcm[29], mcm[30], mcm[31]), file=fileout)

    Usun = mcm[32] * 1000.
    Vsun = mcm[33] * 1000.
    Wsun = mcm[34] * 1000.
    R0 = mcm[-6]
    V0 = mcm[-5]

    M200 = 4. * np.pi * mcm[26] * mcm[30]**3. * (np.log(1. + mcm[-9] / mcm[30]) - mcm[-9] / (mcm[-9] + mcm[30])) / (1.e10)

    c200 = mcm[-9] / mcm[30]
    rs = mcm[30]

    M_NFW = M200
    rs_NFW = rs
    c_NFW = c200

    fileout.close()

    return M_NFW, rs_NFW, c_NFW
