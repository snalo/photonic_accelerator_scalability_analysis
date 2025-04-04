# =============================================
# Silicon Photonic Accelerator Scaling Simulation
# =============================================
# This script evaluates the scalability of MRR-based optical accelerators.
# Key assumption: This simulation is for DIGITAL accelerators (e.g. PIDCOP),
# where photodetector resolution is effectively 1-bit.
# =============================================

import numpy as np
import pandas as pd

def findPDSensitivity(DR):
    # -----------------------------
    # Calculate PD sensitivity threshold (in dBm) needed for 1-bit detection
    # This is based on Eq. (8) from the referenced paper
    # -----------------------------
    no_of_bits = 1  # Fixed to 1 for digital accelerators
    # For digital accelerators (e.g. PIDCOP), the photodetector only needs to distinguish
    # between two levels (presence/absence of light), so effective resolution is always 1-bit,
    # regardless of overall system bit precision.

    e = 1.6e-19
    KT = 0.0259 * e
    R = 1.2
    Id = 35e-9
    RL = 50
    DR_Hz = DR * 1e9
    RIN = 10**(-140 / 10)

    Pd_range = np.arange(-35, 0, 0.1)
    error_list = []

    for pd_dbm in Pd_range:
        Pd = 10**((pd_dbm - 30) / 10)
        A = R * Pd
        B = 2 * e * (R * Pd + Id)
        C = (4 * KT) / RL
        D = (R**2) * (Pd**2) * RIN
        E = 2 * e * Id + C
        F = DR_Hz / np.sqrt(2)
        no_of_bits_hat = (1 / 6.02) * (20 * np.log10(A / ((np.sqrt(B + C + D) + np.sqrt(E)) * np.sqrt(F))) - 1.76)
        error_list.append(abs(no_of_bits_hat - no_of_bits))  # Match predicted bit resolution with 1b

    best_idx = np.argmin(error_list)
    return Pd_range[best_idx]

def findOptimalParallelism(PLaser, Pd_thresh_dbm, tpe, p_penalty):
    # -----------------------------
    # Evaluate max matrix scale (N or M) before optical output drops below PD sensitivity
    # Based on Equation (11) from the paper
    # -----------------------------
    losses = {
        "Psmf_att": 0,
        "Pec_il": 1.6,
        "Psi_att": 0.3,
        "Pmrm_ip_il": 4,
        "Pmrm_ip_obl": 0.01,
        "Psplitter_il": 0.01,
        "Pmrr_w_il": 0.01,
        "Pmrr_w_obl": 0.01,
        "Pmrr_fltr_il": 0.01,
        "dMRR": 0.02
    }

    if tpe == "N=M":
        maxN = 1
        for n in range(1, 300):
            m = n
            Pout = PLaser - losses["Psmf_att"] - losses["Pec_il"] - (losses["Psi_att"] * n * losses["dMRR"]) - \
                   losses["Pmrm_ip_il"] - (n - 1) * losses["Pmrm_ip_obl"] - \
                   (10 * np.log10(n) + losses["Psplitter_il"] * np.log2(m)) - losses["Pmrr_w_il"] - \
                   (n - 1) * losses["Pmrr_w_obl"] - p_penalty
            if Pout < Pd_thresh_dbm:
                break
            maxN = n
        return maxN, Pout
    elif tpe == "NnotM":
        maxM = 1
        n = 8  #Hardcoded for number of wavelength per waveguide. - corresponds to the bits or AND operations done in parallel
        for m in range(1, 500000):
            Pout = PLaser - losses["Psmf_att"] - losses["Pec_il"] - losses["Psplitter_il"] * np.log2(m) - \
                   10 * np.log10(n) - (losses["Psi_att"] * n * m * losses["dMRR"]) - losses["Pmrr_fltr_il"] - \
                   losses["Pmrm_ip_il"] - losses["Pmrr_w_il"] - losses["Pmrr_fltr_il"] * (m + 1) - \
                   (losses["Psi_att"] * n * m * losses["dMRR"]) - p_penalty
            if Pout < Pd_thresh_dbm:
                break
            maxM = m
        return maxM, Pout

def run_sweep(DR_vals=[1, 5, 10], PLaser=10):
    # -----------------------------
    # Main function to sweep across Tiers and data rates
    # Result: Network scale limit (Max N or M) for each configuration
    # -----------------------------
    Tier_values = {
        'Tier_1': {'org_mode': 'AMW', 'tpe': 'N=M', 'p_penalty': 5.8},
        'Tier_2': {'org_mode': 'AMW', 'tpe': 'NnotM', 'p_penalty': 5.8},
        'Tier_3': {'org_mode': 'MAW', 'tpe': 'N=M', 'p_penalty': 4.8},
        'Tier_4': {'org_mode': 'MAW', 'tpe': 'NnotM', 'p_penalty': 4.8},
        'Tier_5': {'org_mode': 'MWA', 'tpe': 'N=M', 'p_penalty': 1.8},
        'Tier_6': {'org_mode': 'MWA', 'tpe': 'NnotM', 'p_penalty': 1.8}
    }

    results = []
    for tier, props in Tier_values.items():
        for DR in DR_vals:
            pd_thresh = findPDSensitivity(DR)
            maxP, Pout = findOptimalParallelism(PLaser, pd_thresh, props["tpe"], props["p_penalty"])
            results.append({
                "Tier": tier,
                "Topology": props["tpe"],
                "Org Mode": props["org_mode"],
                "Bit Precision (B)": 1,  # Key assumption: 1-bit resolution at PD
                "Data Rate (GS/s)": DR,
                "PD Sensitivity (dBm)": round(pd_thresh, 2),
                "Max N or M": maxP,
                "P_out (dBm)": round(Pout, 2)
            })
    return pd.DataFrame(results)

df = run_sweep()
df.to_csv("results/digital_scal_result_fixed_1bit.csv", index=False)
print("Saved as digital_scal_1bit_results.csv")
