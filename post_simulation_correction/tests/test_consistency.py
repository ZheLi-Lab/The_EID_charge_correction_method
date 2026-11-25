"""Synthetic consistency tests; no AMBER/APBS execution or molecular validation."""
import contextlib
import csv
import importlib.util
import io
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

MODULE_PATH = Path(__file__).resolve().parents[1] / "charge_correction.py"
SPEC = importlib.util.spec_from_file_location("preserved_charge_correction", MODULE_PATH)
cc = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = cc
SPEC.loader.exec_module(cc)


class CorrectionConsistencyTests(unittest.TestCase):
    def setUp(self):
        # Avoid __init__: it invokes ambpdb and writes molecular input files.
        self.calc = cc.ChargeCorrection.__new__(cc.ChargeCorrection)
        self.calc.phys_const = cc.PhysicalConstants()
        self.calc.chg_mol = -1
        self.calc.chg_com = 3
        self.calc.chg_pro = 4
        self.rips = {
            "protein_het": {"RIP_kj": 2.0},
            "ligand_het": {"RIP_kj": -3.0},
            "ligand_hom": {"RIP_kj": -8.0},
        }

    def side(self, is_complex):
        with contextlib.redirect_stdout(io.StringIO()):
            return self.calc.one_side_charge_correction(
                self.rips, 60.0, 5000, is_complex=is_complex
            )

    def test_total_is_sum_of_five_terms(self):
        for is_complex in (True, False):
            result = self.side(is_complex)
            self.assertEqual(len(result), 6)
            self.assertAlmostEqual(result["ddG_total"], sum(
                result[key] for key in result if key != "ddG_total"
            ))

    def test_net_usv_combination_matches_screened_net(self):
        const = self.calc.phys_const
        result = self.side(True)
        delta_q2 = 4**2 - 3**2
        expected = -(const.ELS * const.coulomb_factor / 2) * delta_q2 / (6 * const.eps_het)
        self.assertAlmostEqual(result["ddG_NET"] + result["ddG_USV"], expected)

    def test_binding_is_ligand_minus_complex(self):
        self.calc.complex_box_parms = cc.BoxParameters("60 60 60", 60, 5000)
        self.calc.ligand_box_parms = cc.BoxParameters("40 40 40", 40, 2000)
        with patch.object(self.calc, "cal_complex_rip", return_value=self.rips), \
             patch.object(self.calc, "cal_ligand_rip", return_value=self.rips), \
             contextlib.redirect_stdout(io.StringIO()):
            results = self.calc.calculate_binding_correction()
        for key, value in results["binding"].items():
            self.assertAlmostEqual(value, results["ligand"][key] - results["complex"][key])

    def test_summary_csv_units(self):
        results = {side: self.side(side == "complex")
                   for side in ("complex", "ligand", "binding")}
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "summary.csv"
            with contextlib.redirect_stdout(io.StringIO()):
                self.calc.save_correction_results(results, str(output))
            with output.open(newline="") as handle:
                rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 18)
        self.assertEqual(set(rows[0]), {
            "System", "Correction_Type", "Value_kJ_mol", "Value_kcal_mol"
        })
        for row in rows:
            expected = results[row["System"]][row["Correction_Type"]]
            self.assertAlmostEqual(float(row["Value_kJ_mol"]), expected, delta=0.000051)
            self.assertAlmostEqual(float(row["Value_kcal_mol"]), expected / 4.184, delta=0.000051)


if __name__ == "__main__":
    unittest.main()
