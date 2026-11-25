# Validation status and implementation caveats

## Preservation and scope

The two supplied scientific workflow files are unchanged. Existing tracked files
in the upstream repository are also unchanged. No reference molecular dataset
was supplied with this addition; the existing EID energy CSV examples are not
AMBER topology/coordinate fixtures for this workflow.

## Known limitations requiring review

1. **Box dimensions:** `check_cubic_box()` checks the AMBER periodic box, but
   `get_cubic_box_parms()` estimates dimensions from PDB atom-coordinate spans
   using whitespace-split fields. It then truncates the x span to integer Å for
   single-side corrections. These are not necessarily the periodic box lengths;
   PDB chain-field layouts and coordinate wrapping can also affect this estimate.
2. **Atom selection:** `p.pqr` contains every line starting with `ATOM`, whereas
   `l.pqr` contains every line matching the ligand residue-name substring. If
   `ambpdb` marks the ligand as `ATOM`, it can occur in both selections. Water,
   ions, cofactors, and residue-name substring collisions require inspection.
   Zero-charge rewriting assumes fixed-width PQR fields.
3. **Ligand-side geometry:** only the complex is converted with `ambpdb`.
   The ligand-side APBS calculation reuses ligand coordinates/charges extracted
   from that complex; its box estimate and water count come from the separate
   ligand input. Confirm that this geometry choice is intended.
4. **Water model:** water counting is the number of PDB lines containing `WAT`
   divided by three. The physical constants use dielectric 97 and quadrupole
   trace 0.00764 e nm², documented in the supplied code as TIP3P parameters.
   Other residue names or water models are not automatically recognized.
5. **Parameter coupling:** APBS solvent dielectrics and `PhysicalConstants` are
   separate configuration objects; temperature is also fixed to 300 K in the
   template. Changing one object alone does not synchronize all settings.
6. **DX reader:** assumes fixed header-line positions for grid spacings and
   selects numeric lines by their first character. Other DX layouts, leading
   whitespace, or extra numeric records may be misread.
7. **Neutral ligand:** the effective-radius expression divides by the absolute
   ligand charge, so zero charge is unsupported.
8. **CSV fallback:** `save_details_results()` builds its `data` rows after
   importing optional pandas. Without pandas, the fallback references `data`
   before assignment. For a complete workflow with the preserved script,
   install pandas, even though summary export alone has a working fallback.
9. **Execution:** construction of `ChargeCorrection` immediately writes files
   and invokes `ambpdb`. Runs must be isolated in separate working directories.
   APBS output is captured rather than saved as per-task logs; failed jobs may
   require manual rerunning to inspect diagnostics. CLI `-o` writes an additional
   summary after the default workflow output rather than redirecting all output.
10. **Conventions:** binding is ligand minus complex; all five contributions
    enter `ddG_total`. This is distinct from the EID additional correction.
    Confirm applicability and signs against the exact thermodynamic cycle.

## Suggested scientific acceptance checks before a validated release

- Supply one small, redistributable AMBER complex/ligand fixture and provenance.
- Record Python, NumPy, ParmEd, pandas, AMBER/AmberTools, and APBS versions.
- Verify atom counts, charge sums, radii, solvent/ion exclusion, and no duplicate
  atoms in every generated PQR; check both bound and free ligand geometries.
- Compare coordinate spans and integer lengths against actual periodic-box
  lengths, and assess the effect on each correction term.
- Compare all RIP values and NET/USV/RIP/EMP/DSC terms against an independently
  checked reference, including kJ/kcal conversion and binding sign.
- Check grid-size convergence and available memory before batch calculations.
- Archive the verified inputs, expected outputs, settings, and numerical tolerances.
- Address behavioral/scientific changes in a separate revision with before/after
  results, rather than silently altering this preserved implementation.

The included synthetic unit tests verify only internal consistency, not physical
correctness or agreement with the literature.
