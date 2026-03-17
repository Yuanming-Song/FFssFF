# FFssFF All-Atom MD Simulations

Simulation scripts and analysis code for the paper:

> **Nanoscale Interfacial Organization Governs Maturation and Collapse in Passive Versus Active Condensates**

All-atom molecular dynamics simulations of the self-assembling peptide **FFssFF** (Phe-Phe-Ser-Ser-Phe-Phe), studying aggregation under different pH conditions, solvent environments, and concentrations.

**Code written by Yuanming Song. README and code comments co-written with LLM coding tools.**

## Systems

| System | pH | Environment | # Molecules | Box (Å) | Conc. (mM) | Setup modifications |
|---|---|---|---|---|---|---|
| Acid in water | Acidic | Solvated | 32 | 125 | ~29 | Use `FFssFF_Acid.str`, `FFssFF_only_Acid.pdb` |
| Base in water | Basic | Solvated | 32 | 125 | ~29 | Use `FFssFF_base.str`, `FFssFF_only_base.pdb` |
| Base in vacuum | Basic | Vacuum | 32 | 125 (lattice) | — | Use `nvt_eq.inp`/`nvt01.inp` instead of NPT |
| Base high conc. 1 | Basic | Solvated | 32 | 100 | ~53 | Change `set a 100.` in `setUpMix.tcl` |
| Base high conc. 2 | Basic | Solvated | 32 | 85 | ~87 | Change `set a 85.` in `setUpMix.tcl` |
| Base high conc. 3 | Basic | Solvated | 32 | 74 | ~131 | Change `set a 74.` in `setUpMix.tcl` |

## Prerequisites

- **NAMD 3.0** (CPU for equilibration, GPU for production)
- **VMD 1.9.3+** with `tempoUserVMD`, `tempoUtils`, `solvate`, `autoionize`, `psfgen` packages
- **R** with packages: `bio3d`, `doParallel`, `ggplot2`, `dplyr`, `tidyr`, `cowplot`, `Rcpp`, `sna`, `geometry`, `bigmemory`, `SOMMD`
- **CHARMM36 force field** (`par_all36_cgenff.prm`, `top_all36_cgenff.rtf`, `toppar_water_ions_jaf_cgenff.str`)

## Repository Structure

```
forcefield/         FFssFF-specific force field parameters and input PDBs
simulation/         Template setup, equilibration, production, and SLURM scripts
analysis/           Contact analysis R scripts and output data
  Base_Rscript/     Shared R function library
  output/           Pre-computed .rda data and plots
plotting/           Publication plot generation scripts
```

## Reproducing the Simulations

### 1. System Setup

Edit `simulation/setUpMix.tcl`:
- Set `myPSF` and `myDCD` to point to the single-molecule PDB/DCD (from `forcefield/`)
- Set `topoList` to include the appropriate `.str` file
- Set `set a <box_size>` for the desired concentration
- Set `outname` (e.g., `FFssFF_Base_mix`)

Run:
```bash
vmd -dispdev text -e setUpMix.tcl
```

This generates `{outname}_water_ions.psf`, `.pdb`, and `.cons.pdb`.

### 2. Equilibration

Edit `simulation/npt_eq.inp`:
- Set `Dir` to point to where the setup outputs are
- Set `cellBasisVector` to match box size

Submit (CPU):
```bash
sbatch namd_eq.slurm
```

### 3. Production

Edit `simulation/npt01.inp`:
- Set `Dir` to match

Submit (GPU):
```bash
sbatch namd.cuda.slurm
```

For continuation runs:
```bash
bash conti.sh && sbatch namd.cuda.slurm
```

### 4. Vacuum Simulations

Use `nvt_eq.inp` and `nvt01.inp` instead of the NPT configs. No solvation step — use the lattice PSF/PDB directly.

## Analysis

### Contact Analysis

```bash
# For the 3 original systems (Acid, Base, Dry)
sbatch analysis/run_contact_analysis.slurm

# For the concentration sweep (high1, high2, high3)
sbatch analysis/run_contact_analysis_conc.slurm
```

### Plotting

Open `plotting/plot_contact_analysis.Rmd` in RStudio or run:
```bash
Rscript -e "rmarkdown::render('plotting/plot_contact_analysis.Rmd')"
```

## Simulation Parameters

| Parameter | Value |
|---|---|
| Force field | CHARMM36 + CGenFF (FFssFF custom) |
| Water model | TIP3P |
| Temperature | 297.15 K |
| Pressure | 1 atm (NPT) |
| Timestep | 2 fs (production), 1 fs (equilibration) |
| Electrostatics | PME (order 6, 1.0 Å grid) |
| Cutoff | 12 Å, switching from 10 Å |
| Thermostat | Langevin (γ = 1 ps⁻¹) |
| Barostat | Langevin piston (200 ps period, 50 ps decay) |
| Ionic strength | 0.250 M NaCl |
| Output frequency | 10,000 steps (20 ps) |
