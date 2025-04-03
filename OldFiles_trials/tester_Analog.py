import numpy as np
import pandas as pd
from pathlib import Path

#no_of_bits = [4, 8, 12] #Fix ni/p = 4bits
# Error Note:
#ScalabilityMMA covers the Input-Weight Modulation and the final aggregation
#There will loss reduction as cross-talk

def findPDSensitivity(no_of_bits, DR):
    e = 1.6*10**(-19)  # columbs
    KT = 0.0259*e  # columbs*Volt
    R = 1.2  # A/W
    Id = 35*10**(-9)  # A
    RL = 50  # ohm
    DR = DR*10**(9)  # Bits/s
    RIN = 10**(-140/10)  # power ratio/ Hz

    Pd_range = np.arange(-30, 2, 0.1)
    error_list = []
    pd_list = []
    for pd_dbm in Pd_range:
        Pd = 10**((pd_dbm-30)/10)  # W
        A = R*Pd
        B = 2*e*(R*Pd+Id)
        C = (4*KT)/RL
        D = (R**2)*(Pd**2)*RIN
        E = 2*e*Id + C
        F = DR/np.sqrt(2)
        no_of_bits_hat = (
            1/6.02)*(20*np.log10(A/((np.sqrt(B+C+D)+np.sqrt(E))*np.sqrt(F)))-1.76)
        error = abs(no_of_bits_hat - no_of_bits)
        error_list.append(error)
        pd_list.append(pd_dbm)
    min_error = min(error_list)
    min_error_idx = error_list.index(min_error)
    pd = pd_list[min_error_idx]
    print(f"Calculated PD Sensitivity for {DR/10**9}hz is {pd}  ")
    return pd  # PD_OPT

def findN_M_PLaser(pd, N):
    error_list = []
    N_list = []
    M_list = []
    PLaser_list = []
    Pd_dbm = pd  # dBm
    Pd = 10**((Pd_dbm)/10)
    M_range = np.arange(1, 1600, 1) #No of DPEs
    for M in M_range:
        eta_WG = 10 ** ((0.3 * (0.02) * N) / 10)
        eta_SMF = 1
        eta_EC = 10 ** (1.6 / 10)
        eta_WPE = 0.1
        IL_MRM = 10 ** (4 / 10)
        IL_MRR = 10 ** (0.01 / 10)
        OBL_MRM = 10 ** (0.01 / 10)
        OBL_MRR = 10 ** (0.01 / 10)
        IL_Penalty = 10 ** ((1.8) / 10) #4.8-3dB since limited error due to crosstalk
        EL_splitter = 10 ** ((0.01) / 10)
        d_MRR = 0.02  # mm
        eta_EC = 1 / eta_EC
        IL_MRM = 1 / IL_MRM
        OBL_MRM = 1 / OBL_MRM
        EL_splitter = 1 / EL_splitter
        IL_MRR = 1 / IL_MRR
        OBL_MRR = 1 / OBL_MRR
        IL_Penalty = 1 / IL_Penalty
        A = M * eta_WG
        B = eta_SMF * eta_EC * IL_MRM * \
            (OBL_MRM ** (N - 1)) * (EL_splitter ** (np.log2(M)))
        C = Pd
        D = eta_WPE * IL_MRR * (OBL_MRR ** (N - 1)) * IL_Penalty
        PLaser = (A / B) * (C / D)
        PLaser_dbm = 10*np.log10(PLaser)
        PLaser_dbm = PLaser_dbm/2
        error_list.append(abs(10-PLaser_dbm)) #subtracts 10 - calculated PLaser
        N_list.append(N)
        M_list.append(M)
        PLaser_list.append(PLaser_dbm)
    min_error_idx = error_list.index(min(error_list)) #find the index of the min error list
    print("Minimum Error ", min(error_list))
    N = N_list[min_error_idx]
    M = M_list[min_error_idx]
    PLaser = PLaser_list[min_error_idx]
    print("Calculated N value", N)
    print("Calculated M value", M)
    return N, M, PLaser

def findOpticalRecievedPower(N, M, PLaser):
    # # Optical Power Calculation
    Psmf_att = 0
    Pec_il = 1.6
    Psi_att = 0.3
    Pmrm_ip_il = 4
    Pmrm_ip_obl = 0.01
    Psplitter_il = 0.01
    Pmrr_w_il = 0.01
    Pmrr_w_obl = 0.01
    p_penalty = 1.8 #4.8-3dB since limited error due to crosstalk
    Pout = PLaser - Psmf_att - Pec_il - Psi_att - Pmrm_ip_il - (N-1)*Pmrm_ip_obl - (
        10*np.log10(M)+Psplitter_il*np.log2(M)) - Pmrr_w_il - (N-1)*Pmrr_w_obl - p_penalty
    return Pout

bits_range = [4, 8, 12 ] # May be removed
DR_range = [1, 5, 10, 20]
N_range = [1, 4, 9, 16]
result_list = []
for no_of_bits in bits_range:
    for DR in DR_range:
        for N in N_range:
            result = {}
            Pd = findPDSensitivity(no_of_bits, DR)
            N, M, PLaser = findN_M_PLaser(Pd, N)
            Pout = findOpticalRecievedPower(N, M, PLaser)
            result['no_of_bits'] = no_of_bits
            result['DR'] = DR
            result['N'] = N
            result['M'] = M
            result['Pout'] = Pout
            # result['no_of_bits'] = no_of_bits
            print(result)
            result_list.append(result)
df = pd.DataFrame(result_list)
filepath = Path('MMA_tester.xlsx')
filepath.parent.mkdir(parents=True, exist_ok=True)
df.to_excel(filepath)