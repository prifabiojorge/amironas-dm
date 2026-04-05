# Allosteric Inhibition of 5-LOX by Pentacyclic Triterpene Analogs — Molecular Dynamics Analysis

This repository contains the analysis scripts, processed data, and publication figures from the study on allosteric inhibition of human 5-lipoxygenase (5-LOX) by pentacyclic triterpene analogs, evaluated through molecular docking, 100 ns molecular dynamics (MD) simulations, and MM/GBSA binding free energy calculations.

## Contents

- `data/`: Processed simulation data (RMSD, RMSF, Rg, SASA, and MM/GBSA results)
- `scripts/`: Python scripts for data analysis and figure generation
- `figures/`: Publication-quality figures

## Ligands Evaluated

| Compound | Description | ΔG_GB (kcal/mol) |
|----------|-------------|-------------------|
| AKBA (AF7) | Co-crystallized control (PDB: 6NCF) | −33.64 ± 4.58 |
| PubChem CID 101685610 | Optimized derivative (**lead candidate**) | −33.03 ± 3.24 |
| PubChem CID 101663382 | Optimized derivative | −31.75 ± 3.48 |
| α-Amyrone (CID 73170) | Unsubstituted ursane scaffold | −27.82 ± 4.08 |
| 6-Epicyclomusalenone (EPIC) | Natural triterpene (COCONUT) | −29.71 ± 31.72 |

## Requirements

- Python ≥ 3.9 with NumPy, SciPy, Matplotlib
- GROMACS 2025.3 (for reproducing simulations)

## Simulation Details

### Software & Force Fields
- **MD Engine**: GROMACS 2025.3 (Coaraci HPC Cluster)
- **Protein force field**: CHARMM36 all-atom
- **Ligand parameterization**: ACPYPE (GAFF + RESP charges)
- **Water model**: TIP3P
- **Docking**: SwissDock
- **MM/GBSA**: gmx_MMPBSA v1.5.0.3

### Protein Target
- **Structure**: Human stable 5-lipoxygenase (PDB ID: 6NCF)
- **Allosteric site**: C2-like / catalytic domain interface (~30 Å from catalytic iron)
- **Grid coordinates**: x=11.6, y=−23.38, z=−18.01 (30 × 28 × 26 Å)

### MD Protocol
- Energy minimization (steepest descent)
- NVT equilibration (2 stages, V-rescale thermostat, 298.15 K)
- NPT equilibration (2 stages, Parrinello-Rahman barostat, 1.0 bar)
- Production: 100 ns unconstrained NPT ensemble
- Trajectory frames: 1001 (every 100 ps)

### MM/GBSA Parameters
- Frames: 34 representative frames across equilibrated trajectory
- Solvent model: Generalized Born (GB) implicit solvent
- Per-residue decomposition: van der Waals + electrostatic contributions

## Data Organization

### `data/mmpbsa/`
MM/GBSA binding free energies and per-residue decomposition for all 5 systems:
- `*_FINAL.dat` — Total binding free energies (ΔG_GB, ΔE_VDW, ΔE_ELE, ΔG_solv)
- `*_DECOMP.dat` — Per-residue energy decomposition
- `*_FRAME.csv` — Frame-by-frame energies
- `*_DECOMP_FRAME.csv` — Frame-by-frame per-residue decomposition

### `data/rmsd/`, `data/rmsf/`, `data/rg/`, `data/sasa/`
Processed structural analysis data and figures for each system.

## Analysis Scripts

| Script | Purpose |
|--------|---------|
| `scripts/analyze_rg_complete_CORRECTED.py` | Radius of gyration analysis with nm→Å conversion |
| `scripts/analyze_rmsf_complete.py` | Root mean square fluctuation analysis |
| `scripts/analyze_sasa_complete.py` | Solvent accessible surface area analysis |
| `scripts/gen_fig6_decomp.py` | Per-residue energy decomposition figure (Figure 6) |
| `scripts/extract_sim_data.py` | Extract raw simulation data from GROMACS outputs |

## Key Findings

- **PubChem CID 101685610** achieved binding free energy comparable to the co-crystallized AKBA control (−33.03 vs −33.64 kcal/mol)
- Superior dynamic stability with lowest ligand RMSD (1.53 ± 0.43 Å vs 2.72 ± 0.84 Å for AKBA)
- Efficient binding strategy: minimal electrostatic interactions compensated by reduced desolvation penalty
- α-Amyrone and EPIC showed partial/complete dissociation from the allosteric pocket

## License

This work is licensed under a Creative Commons Attribution 4.0 International License (CC-BY 4.0). See the `LICENSE` file for details.

## Citation

If you use this data or code in your research, please cite the original manuscript (currently under preparation) and this repository.

## Contact

For questions regarding this reproducibility package, please contact the corresponding author.

---

*This package was prepared in support of manuscript reproducibility requirements.*
