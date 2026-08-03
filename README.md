# Reproducibility Package

This repository provides the post-processing scripts and reference figures required to reproduce the results presented in the manuscript:

**“Excitation of Plasma Waves by Two-Stream Instability in Unmagnetized Plasma”**

The package is designed to support the computational reproducibility of the published results.

---

## Data Availability

The full Particle-in-Cell (PIC) simulation datasets associated with this study are archived on Zenodo:

**DOI:** https://doi.org/10.5281/zenodo.20202292

The Zenodo archive contains the following simulation datasets:

```text
yaoxpic_v25_counter_1.zip
yaoxpic_v25_counter_2.zip
yaoxpic_v25_counter_3.zip
```

Due to file size constraints, the large simulation datasets are hosted on Zenodo rather than in this GitHub repository.

---

## Repository Contents

This GitHub repository contains:

- Python scripts used to generate all figures in the manuscript
- Auxiliary Python modules used for numerical analysis and wave calculations
- Pre-generated figures for reference and comparison

### Repository Structure

```text
scripts/   → Python scripts for figure generation
fig/       → Pre-generated figures (PDF/PNG)
data/      → Local directory for downloaded simulation datasets
```

The three simulation datasets correspond to:

```text
data/yaoxpic_v25_counter_1/
data/yaoxpic_v25_counter_2/
data/yaoxpic_v25_counter_3/
```

Each dataset contains:

- Particle distribution data
- Time-series diagnostics
- Spectral analysis outputs

All simulation outputs are stored in **HDF5 (`.h5`) format**.

---

## Setup Instructions

### Step 1. Clone the repository

```bash
git clone https://github.com/yaoxin2026/Two-Stream-Instability-Study.git
cd Two-Stream-Instability-Study
```

---

### Step 2. Install required Python packages

```bash
pip install -r requirements.txt
```

Required packages:

```text
numpy
scipy
matplotlib
h5py
```

---

### Step 3. Download simulation datasets from Zenodo

Download the following files from the Zenodo archive:

```text
yaoxpic_v25_counter_1.zip
yaoxpic_v25_counter_2.zip
yaoxpic_v25_counter_3.zip
```

Create a local data directory:

```bash
mkdir -p data
```

Extract the datasets into the `data/` directory:

```bash
unzip yaoxpic_v25_counter_1.zip -d ./data
unzip yaoxpic_v25_counter_2.zip -d ./data
unzip yaoxpic_v25_counter_3.zip -d ./data
```

After extraction, the expected directory structure is:

```text
Two_Stream_Instability_Study/
├── data/
│   ├── yaoxpic_v25_counter_1/
│   ├── yaoxpic_v25_counter_2/
│   └── yaoxpic_v25_counter_3/
├── fig/
├── scripts/
├── README.md
└── requirements.txt
```


---

## Reproducing the Figures


Please generate data about wave branches before any operation.

```bash
cd scripts
python Data_wave_equation_three_runs.py
```



Each figure in the manuscript can be reproduced independently using the corresponding Python script.

Example:

```bash
cd scripts
python Figure01_particle_vspace.py
```


The generated figures can be compared directly with the reference figures provided in the `fig/` directory.

---

## Helper Modules

The repository includes auxiliary Python modules:

```text
YaoxPy_PIC.py
YaoxPy_Wave_Equations.py
YaoxPy_Import_CWD.py
```

These modules provide shared numerical routines, wave-equation solvers, and utility functions used by the figure-generation scripts.

---

## Reproducibility Notes

- No additional preprocessing is required once the simulation datasets are downloaded.
- The workflow is deterministic given identical inputs.
- All scripts read simulation outputs directly from HDF5 files.
- Pre-generated figures are included for validation and comparison.

---

## Contact

For questions regarding the simulation datasets or scripts, please contact:

**Xin Yao**  
Email: yaoxin@nssc.ac.cn

---

## Citation

If you use this repository or the associated simulation datasets, please cite:

```text
Yao, X. (2026).
Reproducibility Package for “Excitation of Plasma Waves by Two-Stream Instability in Unmagnetized Plasma”.
Zenodo. https://doi.org/10.5281/zenodo.20202292
```