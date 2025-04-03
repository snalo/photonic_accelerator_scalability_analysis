import numpy as np


def findPDSensitivity(no_of_bits, DR):
    e = 1.6 * 10 ** (-19)  # columbs
    KT = 0.0259 * e  # columbs*Volt
    R = 1.2  # A/W
    Id = 35 * 10 ** (-9)  # A
    RL = 50  # ohm
    DR = DR * 10 ** (9)  # Bits/s
    RIN = 10 ** (-140 / 10)  # power ratio/ Hz

    Pd_range = np.arange(-35, 0, 0.1)
    error_list = []
    pd_list = []
    for pd_dbm in Pd_range:
        Pd = 10 ** ((pd_dbm - 30) / 10)  # W
        A = R * Pd
        B = 2 * e * (R * Pd + Id)
        C = (4 * KT) / RL
        D = (R ** 2) * (Pd ** 2) * RIN
        E = 2 * e * Id + C
        F = DR / np.sqrt(2)
        no_of_bits_hat = (
                                 1 / 6.02) * (
                                 20 * np.log10(A / ((np.sqrt(B + C + D) + np.sqrt(E)) * np.sqrt(F))) - 1.76)
        error = abs(no_of_bits_hat - no_of_bits)
        error_list.append(error)
        pd_list.append(pd_dbm)

    min_error = min(error_list)
    min_error_idx = error_list.index(min_error)
    pd_dbm = pd_list[min_error_idx]
    print("*******Calculated PD Sensitivity*****", pd_dbm)
    return pd_dbm


import numpy as np


def findOptimalN(PLaser, pd_dbm, arch=""):
    # # Optical Power Calculation
    Psmf_att = 0
    Pec_il = 1.6
    Psi_att = 0.3
    Pmrm_ip_il = 4
    Pmrm_ip_obl = 0.01
    Psplitter_il = 0.01
    Pmrr_w_il = 0.01
    Pmrr_fltr_il = 0.01
    Pmrr_w_obl = 0.01
    dMRR = 0.02;
    if arch == 'MAW':
        p_penalty = 4.8  #
    elif arch == 'AMW':
        p_penalty = 5.8
    elif arch == 'MWA':
        p_penalty = 1.8
    # n_range = 300
    # N = range(1,n_range)
    # Pout =  0
    # maxN =1
    # for n in N:
    #     m = n
    #     Pout = PLaser - Psmf_att - Pec_il - (Psi_att*n*dMRR) - Pmrm_ip_il - (n-1)*Pmrm_ip_obl - (10*np.log10(n)+Psplitter_il*np.log2(m)) -Pmrr_w_il -(n-1)*Pmrr_w_obl - p_penalty
    #     if Pout<pd_dbm:
    #       break
    #     else:
    #       maxN = n
    # for M not equal to N
    m_range = 500000
    M = range(1, m_range)
    Pout = 0
    maxM = 1
    for m in M:
        n = 4
        # Pout = PLaser - Psmf_att - Pec_il - Psplitter_il * np.log2(m) - 10 * np.log10(n) - \
        #        (Psi_att * n * m * dMRR) - Pmrr_fltr_il - Pmrm_ip_il - \
        #        Pmrr_w_il - Pmrr_fltr_il - Pmrr_fltr_il * (m - 1) - Psi_att * n * m * dMRR - p_penalty
        #Modification from Dr. TK
        Pout = PLaser - Psmf_att - Pec_il - Psplitter_il * np.log2(m) - 10 * np.log10(n) - (
                    Psi_att * n * m * dMRR) - Pmrr_fltr_il - Pmrm_ip_il - Pmrr_w_il - Pmrr_fltr_il - Pmrr_fltr_il * (
                           m - 1) - (Psi_att * n * m * dMRR) - p_penalty
        if Pout < pd_dbm:
            break
        else:
            maxM = m
    print("*******Calculated Max Supported N*****", maxM)
    return maxM, Pout


import pandas as pd

bits_range = [4]
DR_range = [1, 5, 10]
PLaser = 10  # dB
arch = 'MWA'
result_list = []
for no_of_bits in bits_range:
    for DR in DR_range:
        result = {}
        Pd_dbm = findPDSensitivity(no_of_bits, DR)
        M, Pout = findOptimalN(PLaser, Pd_dbm, arch)
        # result['PD_Sensitivity'] = Pd
        # result['N'] = N
        result['M'] = M
        # result['PLaser'] = PLaser
        result['Pout'] = Pout
        result['no_of_bits'] = no_of_bits
        result['DR'] = DR
        result_list.append(result)
df = pd.DataFrame(result_list)
df.to_csv('MWA_N_Recieved_Power.csv')