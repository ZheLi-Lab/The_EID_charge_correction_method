# Post-simulation analytical charge correction

[中文说明](README.zh-CN.md) · [Implementation notes and validation limits](VALIDATION.md)

This directory supplements the original repository with the supplied
post-simulation analytical correction workflow. **It is separate from the EID
fitting and `eid_additional_correction.py` workflows.** Do not automatically add
both corrections to the same free-energy estimate: the appropriate correction
and sign depend on the thermodynamic cycle and the underlying calculation.

The supplied `charge_correction.py` and `apbs_tmp.in` are preserved byte-for-byte.
The organization and documentation do not constitute scientific validation.

## Contents

- `charge_correction.py`: AMBER input preparation, APBS calculations, residual
  integrated potentials (RIP), and NET/USV/RIP/EMP/DSC correction terms.
- `apbs_tmp.in`: APBS template; must remain beside the Python script.
- `requirements.txt`: Python packages for the complete workflow (not a version-locked environment).
- `tests/`: synthetic consistency checks without AMBER/APBS or molecular data.
- `VALIDATION.md`: important assumptions, implementation limitations, and release checks.

## Dependencies

Use Python 3.9+ as a recommended environment (older versions have not been tested).

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r post_simulation_correction/requirements.txt
```

The dependency list includes pandas: although summary export has a plain-CSV
fallback, the preserved detailed-export fallback is broken (see `VALIDATION.md`).

Install these external executables separately and make them available on `PATH`:

- `ambpdb` from AMBER/AmberTools, for PQR conversion;
- `apbs`, for Poisson–Boltzmann potential grids.

Check availability with `command -v ambpdb` and `command -v apbs`.
Five APBS tasks run with up to five worker processes by default, each using a
257 × 257 × 257 grid. Plan memory accordingly. The CLI has no worker-count option.

## Inputs and run example

Prepare matched solvated AMBER topology/coordinate pairs for the complex and
ligand-in-solvent systems. Both boxes must be approximately cubic. Coordinates
must be readable by ParmEd; the complex coordinates must also be a text format
accepted by `ambpdb` through stdin (a binary trajectory is not suitable).

From the repository root, use a **new, empty working directory for each run**:

```sh
mkdir -p runs/system_01
cd runs/system_01
# Replace /absolute/path/to/inputs with your own data directory.
python ../../post_simulation_correction/charge_correction.py \
  /absolute/path/to/inputs/complex.rst7 \
  /absolute/path/to/inputs/complex.prmtop \
  /absolute/path/to/inputs/ligand.rst7 \
  /absolute/path/to/inputs/ligand.prmtop \
  -1 3 -r MOL -o results.csv
```

Argument order: `complex_crd complex_top ligand_crd ligand_top chg_mol chg_com`.
Here `-1` and `3` are illustrative ligand and **solute complex** charges, in
units of elementary charge. They are not the net charges of neutralized solvated
boxes. The protein charge used internally is `chg_com - chg_mol`.
`-r` identifies the ligand residue name (default `MOL`). A neutral ligand is not
supported by the effective-radius expression. Run `--help` for CLI syntax.

**Before a production calculation**, read `VALIDATION.md` and inspect the generated
PQR files to ensure that protein/ligand atoms are selected exactly once and that
solvent and ions are excluded as intended. The current selection does not
ensure this automatically.

## Results and conventions

All outputs are written under the current working directory and may overwrite
existing files:

| Path | Content |
| --- | --- |
| `complex.pdb`, `ligand.pdb` | Structures used for box-span and water-count estimates |
| `com.pqr`, `p.pqr`, `l.pqr`, `p0.pqr`, `l0.pqr` | Intermediate charge/radius files |
| `complex/` | Three APBS inputs, potential grids, PQR files, and RIP summaries |
| `ligand/` | Two APBS inputs, potential grids, PQR files, and RIP summaries |
| `charge_corrections.csv` | Default summary for complex, ligand, and binding |
| `charge_corrections_details.csv` | Per-side RIP values, box lengths, charges, radii, and correction terms |
| `results.csv` (if `-o results.csv`) | Additional summary copy; `-o` does not suppress the default files |

The summary contains `System`, `Correction_Type`, `Value_kJ_mol`, and
`Value_kcal_mol`. Terms are `ddG_NET`, `ddG_USV`, `ddG_RIP`, `ddG_EMP`,
`ddG_DSC`, and their sum `ddG_total`.
The implementation defines every binding-level term as **ligand side minus
complex side**. Verify that this matches your annihilation/decoupling and
binding-free-energy sign conventions before applying it.

## Tests

From the repository root:

```sh
python -m unittest discover -s post_simulation_correction/tests -v
```

These tests exercise algebraic consistency and CSV output using synthetic RIP
values. They do not validate PQR selection, APBS execution, or agreement with
published/reference molecular results.

## Scientific context

The terms follow the analytical finite-size charge-correction framework discussed
by Rocklin et al., *Calculating the binding free energies of charged species based
on explicit-solvent simulations employing lattice-sum methods: An accurate
correction scheme for electrostatic finite-size effects*, J. Chem. Phys. **139**,
184103 (2013), https://doi.org/10.1063/1.4826261.
This citation provides context, not a claim that this implementation has been
independently reproduced against that paper. The original repository's EID
workflow and this directory should be described separately in publications.
