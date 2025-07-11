# The_EID_charge_correction_method
The Electrostatic Interaction Decoupling (EID) method is a novel charge correction approach for alchemical free energy calculations of charged ligand-protein systems.
## Installation
```sh
git clone https://github.com/ZheLi-Lab/The_EID_charge_correction_method.git
cd The_EID_charge_correction_method
pip install -e .
```
Dependencies
- Python 3.7+
- NumPy
- Pandas
- SciPy
- Matplotlib (for plotting)
- pymbar (for BAR calculations)
- OpenMM (for simulation data generation)

## Usage
1. Step 1: Generated the Simulation Energy Output Files

Collect the potential energy data from your alchemical simulations. These should be organized by λ windows and contain the electrostatic energy values for each frame.

2. Step 2: Run the EID Analysis
```sh
python EID_charge_correction_analysis.py
```
This script performs multi-window joint fitting of the electrostatic energy data to extract the quadratic coefficients (a, b, c) for each frame, and then uses the linear coefficient (b) to recalculate corrected energy values.

3. Step 3: Analyze the Corrected Free Energy Results
```sh
python AlchemConvTools/one_end_fe_aly.py -i input_aly.txt
```
## Example input_aly.txt Configuration:
```sh
###Input file format:
[Basic_settings]
simulation_software = openmm
file_directory = "."
file_prefix = 'state_ΔUs'
file_suffix = '.csv'
subsample = False
fraction = 1.0
output_csv_filename = 'free_ene.csv'
energy_unit = 'kcal/mol'
calculation_windows = all
std_mode = bar_std
plot_du = False
```
