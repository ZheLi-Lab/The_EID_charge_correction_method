# The_EID_charge_correction_method

This repository contains the implementation of the Electrostatic Interaction
Decoupling (EID) workflow for charge-changing alchemical free energy
calculations.

In the current EID formulation, the strict statement applies to the isolation
of the linear term from the quadratic fit, whereas the final electrostatic free
energy is obtained by adding a residual offset-removal analytical correction to
the linearly separated periodic-boundary-condition (PBC) contribution. The
practical workflow therefore contains two parts:

1. extract the linearly separated PBC contribution from the quadratic fit to
   the electrostatic energy;
2. add the EID additional analytical correction term, which removes the
   residual PBC-induced potential-offset contribution.

This repository therefore provides both the fitting workflow and a minimal
script for the additional EID correction term.

## Installation

```sh
git clone https://github.com/ZheLi-Lab/The_EID_charge_correction_method.git
cd The_EID_charge_correction_method
```

Dependencies for the fitting workflow:

- Python 3.7+
- NumPy
- Pandas
- SciPy
- Matplotlib
- pymbar

## Usage

### 1. Prepare electrostatic energy data

Collect the potential-energy data from the electrostatic decoupling windows.
The example files in `example/` show the expected file layout.

### 2. Perform multi-window quadratic fitting and extract the linear term

```sh
python EID_charge_correction_analysis.py
```

This script fits the electrostatic energy to `U = A * lambda^2 + B * lambda + C`
for each frame and writes the linearly separated PBC dataset derived from the
fitted `B` term.

### 3. Compute the linearly separated PBC free energy

```sh
python AlchemConvTools/one_end_fe_aly.py -i example/input_aly.txt
```

This step applies BAR to the linearly separated dataset and yields the
linearly separated PBC contribution.

### 4. Add the residual offset-removal analytical correction

```sh
python eid_additional_correction.py \
  --ligand-charge -1 \
  --environment-charge 0 \
  --box-length 60 \
  --system-name ligand_in_complex
```

The script writes a one-row CSV file named
`eid_additional_correction.csv` by default and prints the correction in both
kJ/mol and kcal/mol. In the standard EID setup, `environment_charge` should be
zero. The final EID free energy is obtained by adding this term to the
linearly separated PBC contribution from step 3.

If you want the binding-level additional correction directly, you can provide
the complex-side and ligand-side box lengths together:

```sh
python eid_additional_correction.py \
  --ligand-charge -1 \
  --complex-box-length 60 \
  --ligand-box-length 40 \
  --system-name example_binding
```

In this mode, the output CSV contains the complex-side term, the ligand-side
term, and the final binding-side difference for the additional correction.

## Example `input_aly.txt`

```ini
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
