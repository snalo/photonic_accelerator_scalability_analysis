# Photonic Accelerator Scalability Analysis

This project provides a detailed simulation framework to analyze the **scaling limits of silicon photonic-based accelerators**, specifically those using **Microring Resonator (MRR)** and **wavelength-division multiplexing (WDM)** architectures.

The goal is to determine how large an optical matrix (in terms of `N × N` or `N × M`) can scale before the **optical power at the photodetector (PD)** drops below the **sensitivity threshold** required for reliable signal detection.

---

## Background

This work builds upon the key theoretical models and assumptions from the following paper:

> **Scaling up silicon photonic-based accelerators: Challenges and opportunities**  
> **Authors:** M. A. Al-Qadasi, L. Chrostowski, B. J. Shastri, and S. Shekhar  
> **Published in:** APL Photonics, Volume 7, Number 2, February 2022  
> [DOI: 10.1063/5.0070992](http://dx.doi.org/10.1063/5.0070992)

### BibTeX
```bibtex
@article{Al_Qadasi_2022,
   title={Scaling up silicon photonic-based accelerators: Challenges and opportunities},
   volume={7},
   ISSN={2378-0967},
   url={http://dx.doi.org/10.1063/5.0070992},
   DOI={10.1063/5.0070992},
   number={2},
   journal={APL Photonics},
   publisher={AIP Publishing},
   author={Al-Qadasi, M. A. and Chrostowski, L. and Shastri, B. J. and Shekhar, S.},
   year={2022},
   month=feb
}
```

---

## Repository Structure

```
📁 OldFiles_trials                     → Archived or preliminary versions  
📁 results                            → Output CSVs and plots from simulations  
📄 Final_Analog_Scalability.xlsx      → Manually processed analog results  
📄 digital_scal_1bit.py               → Digital-only version (1-bit PD resolution)  
📄 new_digital_dpu_scal.py            → Optimized digital DPU scaling version  
📄 old_Digital_Arch_Scalability.py    → Legacy logic  
📄 scalability_fullSweep_multi_bit.py → Analog + Digital sweep (1–4 bits)  
📄 README.md                          → You’re reading it!  
```

---

## Simulation Workflow

Each script simulates the scalability of photonic accelerators under different **bit resolutions**, **data rates**, and **architectural tiers** (e.g., AMW, MAW, MWA). The simulation:

1. **Calculates the required PD sensitivity** using the signal-to-noise ratio (SNR) model.  
2. **Simulates the optical link budget** using realistic losses across waveguides, microrings, modulators, and splitters.  
3. Determines the **maximum scalable matrix size** (`N` or `M`) before the output power dips below sensitivity.  

---

## Core Equations Used

### Equation (8): Photodetector Resolution Estimation

Used to compute minimum PD power (sensitivity) required to support a target bit resolution `n`.

![Equation (8): Photodetector Resolution Estimation](/assets/Eqn_8.jpg)

Where:  
- \( R \): PD responsivity  
- \( I_d \): Dark current  
- \( R_L \): Load resistance  
- \( q \): Electron charge  
- \( kT \): Thermal noise factor  
- \( DR \): Data rate in Hz  
- RIN: Relative intensity noise  

---

### Equation (11): Optical Link Budget (MRR-based Implementation)

![Equation (11): Optical Link Budget (MRR-based Implementation)](/assets/Eqn_11.jpg)

This is modeled in the code using a detailed breakdown of losses from:  
- Edge coupling  
- Waveguide attenuation  
- Microring insertion losses  
- WDM splitters and combiners  
- Fabrication detuning penalties  

---

## Output CSVs

Each script generates `.csv` outputs containing:

- **Tier** – Tier label (e.g., Tier_1 to Tier_6)  
- **Topology** – Whether N=M (matrix) or N≠M (broadcast)  
- **Org Mode** – Architecture org mode (AMW, MAW, MWA)  
- **Bit Precision (B)** – 1–4 bits  
- **Accelerator Type** – `"Digital"` if B=1, otherwise `"Analog"`  
- **Data Rate (GS/s)** – Data rates from 1 to 10 GS/s  
- **PD Sensitivity (dBm)** – Required optical power at the photodetector  
- **Max N or M** – Maximum size before crossing power threshold  
- **P_out (dBm)** – Output optical power at that limit  

---

## Example Use Cases

- Evaluate **energy scaling trade-offs** with bit precision and modulation  
- Explore design feasibility across **different photonic org modes**  
- Compare scalability between **analog and digital photonic accelerators**  

---

## How to Run

```bash
python digital_scal_1bit.py               # Run 1-bit digital analysis
python scalability_fullSweep_multi_bit.py # Sweep 1-4 bits for analog/digital
```

Output `.csv` files will be found in the `results/` folder or the same directory as the script.

---

For questions, suggestions, or collaborations, feel free to reach out!
