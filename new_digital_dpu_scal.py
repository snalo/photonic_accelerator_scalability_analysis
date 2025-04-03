
import numpy as np
import pandas as pd

# Digital PIDCOP Model:
# - Fully digital pulse-based operation
# - ADC per channel (high power usage)
# - No DACs used

# Parameters
bits_range = [4]
DR_range = [1, 5, 10]
Tier_values = {
    'Tier_1': {'org_mode': 'AMW', 'tpe': 'N=M'},
    'Tier_2': {'org_mode': 'AMW', 'tpe': 'NnotM'},
    'Tier_3': {'org_mode': 'MAW', 'tpe': 'N=M'},
    'Tier_4': {'org_mode': 'MAW', 'tpe': 'NnotM'},
    'Tier_5': {'org_mode': 'MWA', 'tpe': 'N=M'},
    'Tier_6': {'org_mode': 'MWA', 'tpe': 'NnotM'}
}

ADC_power_per_bit = 1.0  # mW per bit per GS/s
power_budget_mW = 1000  # total power budget

results = []

for tier_name, tier_info in Tier_values.items():
    tpe = tier_info['tpe']
    org_mode = tier_info['org_mode']

    for B in bits_range:
        for DR in DR_range:
            adc_power_per_channel = ADC_power_per_bit * B * DR

            if tpe == 'N=M':
                maxN = 1
                for n in range(1, 300):
                    m = n
                    total_adc_power = n * m * adc_power_per_channel
                    if total_adc_power > power_budget_mW:
                        break
                    maxN = n
                results.append({
                    'Tier': tier_name,
                    'Org Mode': org_mode,
                    'Topology': tpe,
                    'Bit Precision (B)': B,
                    'Data Rate (GS/s)': DR,
                    'Max N or M': maxN,
                    'Total ADC Power (mW)': round(maxN**2 * adc_power_per_channel, 2)
                })
            elif tpe == 'NnotM':
                maxM = 1
                n = 8
                for m in range(1, 50000):
                    total_adc_power = n * m * adc_power_per_channel
                    if total_adc_power > power_budget_mW:
                        break
                    maxM = m
                results.append({
                    'Tier': tier_name,
                    'Org Mode': org_mode,
                    'Topology': tpe,
                    'Bit Precision (B)': B,
                    'Data Rate (GS/s)': DR,
                    'Max N or M': maxM,
                    'Total ADC Power (mW)': round(n * maxM * adc_power_per_channel, 2)
                })

df = pd.DataFrame(results)
df.to_csv('results/digital_pidcop_results.csv', index=False)
