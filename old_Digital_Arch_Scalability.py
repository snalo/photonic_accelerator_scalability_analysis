import numpy as np
import pandas as pd


def findPDSensitivity(no_of_bits, DR):  # Function use the no of bits and DR to return the calculated Pd
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
        # Decomposing the noise term Beta into B, C, D, E
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


def findOptimalN(PLaser, pd_dbm, tpe=''):
    # # Optical Power Calculation
    Psmf_att = 0
    Pec_il = 1.6
    Psi_att = 0.3
    Pmrm_ip_il = 4
    Pmrm_ip_obl = 0.01
    Psplitter_il = 0.01
    Pmrr_w_il = 0.01
    Pmrr_w_obl = 0.01
    Pmrr_fltr_il = 0.01
    dMRR = 0.02;
    if tpe == 'N=M':
        # for M = N
        n_range = 300
        N = range(1, n_range)
        Pout = 0
        maxN = 1
        for n in N:
            m = n
            Pout = PLaser - Psmf_att - Pec_il - (Psi_att * n * dMRR) - Pmrm_ip_il - (n - 1) * Pmrm_ip_obl - (
                    10 * np.log10(n) + Psplitter_il * np.log2(m)) - Pmrr_w_il - (n - 1) * Pmrr_w_obl - p_penalty
            if Pout < pd_dbm:
                break
            else:
                maxN = n
        print("*******Calculated Max Supported N*****", maxN)
        return maxN, Pout
    elif tpe == 'NnotM':
        # for M not equal to N
        m_range = 500000
        M = range(1, m_range)
        Pout = 0
        maxM = 1
        for m in M:
            n = 8
            # Pout = PLaser - Psmf_att - Pec_il - (Psi_att * n * dMRR) - Pmrm_ip_il - \
            #        (10 * np.log10(n) + Psplitter_il * np.log2(m)) \
            #        - Pmrr_w_il - (n * m - 1) * Pmrr_w_obl - p_penalty
            Pout = PLaser - Psmf_att - Pec_il - Psplitter_il * np.log2(m) - 10 * np.log10(n) - (
                    Psi_att * n * m * dMRR) - Pmrr_fltr_il - Pmrm_ip_il - Pmrr_w_il - Pmrr_fltr_il - Pmrr_fltr_il * (
                           m - 1) - (Psi_att * n * m * dMRR) - p_penalty
            if Pout < pd_dbm:
                break
            else:
                maxM = m
        print("*******Calculated Max Supported M*****", maxM)
        return maxM, Pout


#PLaser = 10  # dB Laser Power
bits_range = [1] #Bit Precision
DR_range = [1, 5, 10]
PLaser = 10  # dBm
# org_mode = 'MWA'
result_list = []
for no_of_bits in bits_range:
    for DR in DR_range:
        result = {}
        Pd_dbm = findPDSensitivity(no_of_bits, DR)
        # NorM, Pout = findOptimalN(PLaser, Pd_dbm, org_mode)
        # result['PD_Sensitivity'] = Pd
        Tier_values = {
            'Tier_1': {'org_mode': 'AMW', 'tpe': 'N=M'},
            'Tier_2': {'org_mode': 'AMW', 'tpe': 'NnotM'},
            'Tier_3': {'org_mode': 'MAW', 'tpe': 'N=M'},
            'Tier_4': {'org_mode': 'MAW', 'tpe': 'NnotM'},
            'Tier_5': {'org_mode': 'MWA', 'tpe': 'N=M'},
            'Tier_6': {'org_mode': 'MWA', 'tpe': 'NnotM'}
        }
        # Select the Tier
        selected_Tier = 'Tier_1'

        # Get the values for org_mode and tpe based on the selected Tier
        org_mode = Tier_values[selected_Tier]['org_mode']
        tpe = Tier_values[selected_Tier]['tpe']

        # Define the conditions
        Tier_1 = org_mode == 'AMW' and tpe == 'N=M'
        Tier_2 = org_mode == 'AMW' and tpe == 'NnotM'
        Tier_3 = org_mode == 'MAW' and tpe == 'N=M'
        Tier_4 = org_mode == 'MAW' and tpe == 'NnotM'
        Tier_5 = org_mode == 'MWA' and tpe == 'N=M'
        Tier_6 = org_mode == 'MWA' and tpe == 'NnotM'
        # print("Tier_1:", Tier_2)
        # print("org_mode:", org_mode)
        # print("tpe:", tpe)
        # Check the conditions using if statements
        if selected_Tier == 'Tier_1':
            print("The Selected TPC Organisation is", org_mode, "when", tpe)
            p_penalty = 5.8
            NorM, Pout = findOptimalN(PLaser, Pd_dbm, tpe)
        elif selected_Tier == 'Tier_2':
            print("The Selected TPC Organisation is", org_mode, "when", tpe)
            p_penalty = 5.8
            NorM, Pout = findOptimalN(PLaser, Pd_dbm, tpe)
        elif selected_Tier == 'Tier_3':
            print("The Selected TPC Organisation is", org_mode, "when", tpe)
            p_penalty = 4.8
            NorM, Pout = findOptimalN(PLaser, Pd_dbm, tpe)
        elif selected_Tier == 'Tier_4':
            print("The Selected TPC Organisation is", org_mode, "when", tpe)
            p_penalty = 4.8
            NorM, Pout = findOptimalN(PLaser, Pd_dbm, tpe)
        elif selected_Tier == 'Tier_5':
            print("The Selected TPC Organisation is", org_mode, "when", tpe)
            p_penalty = 1.8
            NorM, Pout = findOptimalN(PLaser, Pd_dbm, tpe)
        else:
            print("The Selected TPC Organisation is", org_mode, "when", tpe)
            p_penalty = 1.8
            NorM, Pout = findOptimalN(PLaser, Pd_dbm, tpe)
        result['org_mode'] = org_mode
        result['NM'] = tpe
        result['N size/M Count'] = NorM
        result['PLaser'] = PLaser
        result['Pout'] = Pout
        result['no_of_bits'] = no_of_bits
        result['DR'] = DR
        #result['Number of Wavelength'] = n
        result_list.append(result)
df = pd.DataFrame(result_list)
df.to_csv('NnotM_N&M Results_8bit.csv')

# NOTE:
# AMW:
# Pout = PLaser - Psmf_att - Pec_il - (Pmrr_fltr_il*n) - (Psplitter_il*np.log2(m)) - 10*np.log10(n) -(Psi_att*2*n*dMRR) - Pmrm_ip_il - (Pmrm_ip_obl*(n-1)) -Pmrr_w_il -(Pmrr_w_obl*(n-1)) - p_penalty
#
# MAW:
# Pout = PLaser - Psmf_att - Pec_il - (Pmrr_fltr_il*n) - 10*np.log10(n) - Pmrm_ip_il - (Pmrm_ip_obl*(n-1)) - (Psi_att*n*dMRR) - Psplitter_il*np.log2(m) - Pmrr_w_il - (Pmrr_w_obl*(n-1)) - (Psi_att*n*dMRR) - p_penalty
#
# MWA:
# Pout = PLaser - Psmf_att - Pec_il - Psplitter_il*np.log2(m) - 10*np.log10(n) - (Psi_att*n*m*dMRR) - Pmrr_fltr_il - Pmrm_ip_il - Pmrr_w_il - Pmrr_fltr_il - Pmrr_fltr_il*(m-1) - Psi_att*n*m*dMRR -  p_penalty
# Pmrr_fltr_il = 0.01
