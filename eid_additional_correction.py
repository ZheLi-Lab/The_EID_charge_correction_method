#!/usr/bin/env python3
"""Compute the residual offset-removal analytical correction used in EID.

This script intentionally keeps only the minimal calculation required for the
extra EID correction term.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Union


@dataclass(frozen=True)
class PhysicalConstants:
    """Constants for the EID residual analytical correction."""

    solvent_dielectric: float = 97.0
    coulomb_factor: float = 138.93545585
    lattice_sum_constant: float = -2.837297

    @property
    def prefactor(self) -> float:
        """Return xi_LS / (8 * pi * eps0) in kJ nm / (mol e^2)."""
        return self.lattice_sum_constant * self.coulomb_factor / 2.0


def calculate_eid_additional_correction(
    ligand_charge: float,
    environment_charge: float,
    box_length_angstrom: float,
    constants: PhysicalConstants,
) -> float:
    """Compute the EID residual analytical correction in kJ/mol.

    This term is added to the linearly separated PBC contribution to obtain the
    final EID free energy.
    """

    if box_length_angstrom <= 0:
        raise ValueError("box_length_angstrom must be positive.")

    box_length_nm = box_length_angstrom / 10.0
    coupled_charge = environment_charge + ligand_charge
    delta_q2 = environment_charge**2 - coupled_charge**2
    dielectric_factor = 1.0 - 1.0 / constants.solvent_dielectric
    return constants.prefactor * dielectric_factor * delta_q2 / box_length_nm


def write_rows(
    output_path: Path, rows: Iterable[Dict[str, Union[float, str]]]
) -> None:
    """Write CSV rows summarizing the calculation."""

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "entry",
                "ligand_charge_e",
                "environment_charge_e",
                "box_length_angstrom",
                "eid_additional_correction_kj_per_mol",
                "eid_additional_correction_kcal_per_mol",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate the residual offset-removal analytical correction in the "
            "EID workflow."
        )
    )
    parser.add_argument(
        "--ligand-charge",
        type=float,
        required=True,
        help="Ligand net charge in elementary charge units.",
    )
    parser.add_argument(
        "--environment-charge",
        type=float,
        default=0.0,
        help=(
            "Net charge of all non-ligand atoms (Q_env). The standard EID "
            "workflow uses 0."
        ),
    )
    parser.add_argument(
        "--box-length",
        type=float,
        help="Cubic box length in angstrom for a single-system calculation.",
    )
    parser.add_argument(
        "--solvent-dielectric",
        type=float,
        default=97.0,
        help="Relative dielectric constant of the solvent.",
    )
    parser.add_argument(
        "--system-name",
        default="system",
        help="Name written to the output summary.",
    )
    parser.add_argument(
        "--complex-box-length",
        type=float,
        help="Cubic box length in angstrom for the complex-side correction.",
    )
    parser.add_argument(
        "--ligand-box-length",
        type=float,
        help="Cubic box length in angstrom for the ligand-side correction.",
    )
    parser.add_argument(
        "--complex-environment-charge",
        type=float,
        default=0.0,
        help=(
            "Environment charge for the complex side. The standard EID setup "
            "uses 0."
        ),
    )
    parser.add_argument(
        "--ligand-environment-charge",
        type=float,
        default=0.0,
        help=(
            "Environment charge for the ligand side. The standard EID setup "
            "uses 0."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("eid_additional_correction.csv"),
        help="Output CSV path.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    constants = PhysicalConstants(solvent_dielectric=args.solvent_dielectric)
    binding_mode = (
        args.complex_box_length is not None or args.ligand_box_length is not None
    )

    if binding_mode:
        if args.complex_box_length is None or args.ligand_box_length is None:
            raise ValueError(
                "--complex-box-length and --ligand-box-length must be provided together."
            )

        complex_correction_kj = calculate_eid_additional_correction(
            ligand_charge=args.ligand_charge,
            environment_charge=args.complex_environment_charge,
            box_length_angstrom=args.complex_box_length,
            constants=constants,
        )
        ligand_correction_kj = calculate_eid_additional_correction(
            ligand_charge=args.ligand_charge,
            environment_charge=args.ligand_environment_charge,
            box_length_angstrom=args.ligand_box_length,
            constants=constants,
        )
        binding_correction_kj = ligand_correction_kj - complex_correction_kj

        write_rows(
            output_path=args.output,
            rows=[
                {
                    "entry": f"{args.system_name}_complex",
                    "ligand_charge_e": args.ligand_charge,
                    "environment_charge_e": args.complex_environment_charge,
                    "box_length_angstrom": args.complex_box_length,
                    "eid_additional_correction_kj_per_mol": complex_correction_kj,
                    "eid_additional_correction_kcal_per_mol": complex_correction_kj
                    / 4.184,
                },
                {
                    "entry": f"{args.system_name}_ligand",
                    "ligand_charge_e": args.ligand_charge,
                    "environment_charge_e": args.ligand_environment_charge,
                    "box_length_angstrom": args.ligand_box_length,
                    "eid_additional_correction_kj_per_mol": ligand_correction_kj,
                    "eid_additional_correction_kcal_per_mol": ligand_correction_kj
                    / 4.184,
                },
                {
                    "entry": f"{args.system_name}_binding",
                    "ligand_charge_e": args.ligand_charge,
                    "environment_charge_e": "",
                    "box_length_angstrom": "",
                    "eid_additional_correction_kj_per_mol": binding_correction_kj,
                    "eid_additional_correction_kcal_per_mol": binding_correction_kj
                    / 4.184,
                },
            ],
        )

        print(
            "Residual offset-removal analytical correction "
            f"(complex): {complex_correction_kj:.6f} kJ/mol"
        )
        print(
            "Residual offset-removal analytical correction "
            f"(ligand): {ligand_correction_kj:.6f} kJ/mol"
        )
        print(
            "Residual offset-removal analytical correction "
            f"(binding): {binding_correction_kj:.6f} kJ/mol"
        )
        print(
            "Residual offset-removal analytical correction (binding): "
            f"{binding_correction_kj / 4.184:.6f} kcal/mol"
        )
    else:
        if args.box_length is None:
            raise ValueError(
                "--box-length is required when complex/ligand box lengths are not provided."
            )

        correction_kj = calculate_eid_additional_correction(
            ligand_charge=args.ligand_charge,
            environment_charge=args.environment_charge,
            box_length_angstrom=args.box_length,
            constants=constants,
        )

        write_rows(
            output_path=args.output,
            rows=[
                {
                    "entry": args.system_name,
                    "ligand_charge_e": args.ligand_charge,
                    "environment_charge_e": args.environment_charge,
                    "box_length_angstrom": args.box_length,
                    "eid_additional_correction_kj_per_mol": correction_kj,
                    "eid_additional_correction_kcal_per_mol": correction_kj / 4.184,
                }
            ],
        )

        print(
            "Residual offset-removal analytical correction: "
            f"{correction_kj:.6f} kJ/mol"
        )
        print(
            "Residual offset-removal analytical correction: "
            f"{correction_kj / 4.184:.6f} kcal/mol"
        )

    print(f"Summary written to {args.output}")


if __name__ == "__main__":
    main()
