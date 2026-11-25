# Post-simulation analytical charge correction

## Purpose

Calculate the NET, USV, RIP, EMP, and DSC analytical charge-correction terms
for complex and ligand-in-solvent systems using AMBER inputs and APBS.
This workflow is separate from the EID fitting and EID additional correction.

## Dependencies

Python 3.9+ with:

```sh
pip install numpy parmed pandas
```

Install AmberTools (providing `ambpdb`) and APBS separately, and ensure both
executables are available on `PATH`. Keep `apbs_tmp.in` beside
`charge_correction.py`.

## Run

Prepare matched AMBER topology and coordinate files for the solvated complex
and ligand systems, both with approximately cubic boxes. Use text coordinates
compatible with `ambpdb`, such as AMBER restart files.

From the repository root:

```sh
mkdir -p runs/system_01
cd runs/system_01
python ../../post_simulation_correction/charge_correction.py \
  /path/to/complex.rst7 /path/to/complex.prmtop \
  /path/to/ligand.rst7 /path/to/ligand.prmtop \
  -1 3 -r MOL -o results.csv
```

Replace the paths and charges with your own values. Argument order:

```text
complex_crd complex_top ligand_crd ligand_top chg_mol chg_com
```

- `chg_mol`: ligand charge; must be nonzero.
- `chg_com`: solute complex charge, excluding solvent and counterions.
- `-r`: ligand residue name (default: `MOL`).
- `-o`: additional summary CSV filename (default: `charge_corrections.csv`).

Use a separate working directory for each run to avoid overwriting outputs.
Inspect generated PQR files to confirm correct protein/ligand selection.

## Output format

### Summary: `charge_corrections.csv`

An 18-row table: three systems (`complex`, `ligand`, `binding`), each with
`ddG_NET`, `ddG_USV`, `ddG_RIP`, `ddG_EMP`, `ddG_DSC`, and `ddG_total`.

```csv
System,Correction_Type,Value_kJ_mol,Value_kcal_mol
```

Values are saved to four decimal places. `ddG_total` is the sum of the five
terms. Each binding-level term is **ligand minus complex**; check that this
matches your thermodynamic-cycle sign convention.

### Details: `charge_corrections_details.csv`

Two rows (`complex` and `ligand`) containing water counts, box lengths,
protein/ligand charges, RIP values, effective ligand radii, and individual
correction terms in kJ/mol.

### Other files

- The filename specified by `-o` is an additional summary copy; the default
  summary and details files are still generated.
- `complex/` and `ligand/` contain APBS inputs, DX potential grids, PQR files,
  and RIP reports.
- Intermediate PDB and PQR files are also written in the working directory.
