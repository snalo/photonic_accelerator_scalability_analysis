
import numpy as np
import pandas as pd

def findPDSensitivity(no_of_bits, DR):
    e = 1.6 * 10 ** (-19)
    KT = 0.0259 * e
    R = 1.2
    Id = 35 * 10 ** (-9)
    RL = 50
    DR = DR * 10 ** (9)
    RIN = 10 ** (-140 / 10)

    Pd_range = np.arange(-35, 0, 0.1)
    error_list = []
    pd_list = []
    for pd_dbm in Pd_range:
        Pd = 10 ** ((pd_dbm - 30) / 10)
        A = R * Pd
        B = 2 * e * (R * Pd + Id)
        C = (4 * KT) / RL
        D = (R ** 2) * (Pd ** 2) * RIN
        E = 2 * e * Id + C
        F = DR / np.sqrt(2)
        no_of_bits_hat = (1 / 6.02) * (20 * np.log10(A / ((np.sqrt(B + C + D) + np.sqrt(E)) * np.sqrt(F))) - 1.76)
        error = abs(no_of_bits_hat - no_of_bits)
        error_list.append(error)
        pd_list.append(pd_dbm)

    min_error_idx = np.argmin(error_list)
    return pd_list[min_error_idx]

def findOptimalN(PLaser, pd_dbm, tpe='', p_penalty=1.8):
    Psmf_att = 0
    Pec_il = 1.6
    Psi_att = 0.3
    Pmrm_ip_il = 4
    Pmrm_ip_obl = 0.01
    Psplitter_il = 0.01
    Pmrr_w_il = 0.01
    Pmrr_w_obl = 0.01
    Pmrr_fltr_il = 0.01
    dMRR = 0.02

    if tpe == 'N=M':
        n_range = 300
        maxN = 1
        for n in range(1, n_range):
            m = n
            Pout = PLaser - Psmf_att - Pec_il - (Psi_att * n * dMRR) - Pmrm_ip_il - (n - 1) * Pmrm_ip_obl - \
                   (10 * np.log10(n) + Psplitter_il * np.log2(m)) - Pmrr_w_il - (n - 1) * Pmrr_w_obl - p_penalty
            if Pout < pd_dbm:
                break
            else:
                maxN = n
        return maxN, Pout
    elif tpe == 'NnotM':
        m_range = 500000
        maxM = 1
        for m in range(1, m_range):
            n = 8
            Pout = PLaser - Psmf_att - Pec_il - Psplitter_il * np.log2(m) - 10 * np.log10(n) - \
                   (Psi_att * n * m * dMRR) - Pmrr_fltr_il - Pmrm_ip_il - Pmrr_w_il - Pmrr_fltr_il - \
                   Pmrr_fltr_il * (m - 1) - (Psi_att * n * m * dMRR) - p_penalty
            if Pout < pd_dbm:
                break
            else:
                maxM = m
        return maxM, Pout

# Parameters
bits_range = [4]
DR_range = [1, 5, 10]
PLaser = 10

Tier_values = {
    'Tier_1': {'org_mode': 'AMW', 'tpe': 'N=M', 'p_penalty': 5.8},
    'Tier_2': {'org_mode': 'AMW', 'tpe': 'NnotM', 'p_penalty': 5.8},
    'Tier_3': {'org_mode': 'MAW', 'tpe': 'N=M', 'p_penalty': 4.8},
    'Tier_4': {'org_mode': 'MAW', 'tpe': 'NnotM', 'p_penalty': 4.8},
    'Tier_5': {'org_mode': 'MWA', 'tpe': 'N=M', 'p_penalty': 1.8},
    'Tier_6': {'org_mode': 'MWA', 'tpe': 'NnotM', 'p_penalty': 1.8}
}

results = []

for tier_name, tier_info in Tier_values.items():
    tpe = tier_info['tpe']
    org_mode = tier_info['org_mode']
    p_penalty = tier_info['p_penalty']
    
    for no_of_bits in bits_range:
        for DR in DR_range:
            pd_dbm = findPDSensitivity(no_of_bits, DR)
            NorM, Pout = findOptimalN(PLaser, pd_dbm, tpe, p_penalty)
            results.append({
                'Tier': tier_name,
                'Org Mode': org_mode,
                'Topology': tpe,
                'Bit Precision (B)': no_of_bits,
                'Data Rate (GS/s)': DR,
                'Max N or M': NorM,
                'P_out (dBm)': round(Pout, 2),
                'PD Sensitivity (dBm)': round(pd_dbm, 2)
            })

df = pd.DataFrame(results)
df.to_csv('all_tier_analog_dpu_results.csv', index=False)
