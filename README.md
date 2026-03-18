# FFssFF All-Atom MD Simulations

**Code written by Yuanming Song.**

Simulation scripts and analysis code for the paper:

> **Nanoscale Interfacial Organization Governs Maturation and Collapse in Passive Versus Active Condensates**

## Background

FFssFF (bis(phenylalanyl-phenylalanyl) cystamine) is a short peptide that undergoes pH-dependent phase separation driven by hydrophobic interactions between phenylalanine dipeptides. These droplets can be reversibly formed and dissolved by controlling redox chemistry. All-atom MD simulations were performed to understand the intermolecular interactions within different regions of FFssFF coacervate droplets. Acidic (+2 charged) vs. basic (neutral) conditions were compared to examine whether redox-responsive aggregation behavior is captured by the force field. Basic solvated vs. vacuum simulations serve as proxies for the hydrated outer layer and dense, dehydrated core of FFssFF droplets, as observed in cryo-TEM images. Additional simulations at elevated concentrations probe concentration-dependent aggregation.

## Force Field

FFssFF force field parameters were generated using the CHARMM General Force Field (CGenFF) program version 3.0 with CGenFF version 4.6. The parameter files are provided in `forcefield/`.

## Systems

All production runs performed on NVIDIA A30 GPUs using NAMD 3.0alpha13. All systems contain 32 FFssFF molecules initially arranged on a 5x5x5 lattice.

| System | Charge | Environment | Box (Å) | Conc. (mM) | Speed (ns/day) | Setup modifications |
|---|---|---|---|---|---|---|
| Acid in water | +2 | Solvated | 125 | ~29 | ~20 | Use `FFssFF_Acid.str`, `FFssFF_only_Acid.pdb` |
| Base in water | Neutral | Solvated | 125 | ~29 | ~11 | Use `FFssFF_base.str`, `FFssFF_only_base.pdb` |
| Base in vacuum | Neutral | Vacuum | 125 | — | ~143 | Use `nvt_eq.inp`/`nvt01.inp` instead of NPT |
| Base high conc. 1 | Neutral | Solvated | 100 | ~53 | ~21 | Change `set a 100.` in `setUpMix.tcl` |
| Base high conc. 2 | Neutral | Solvated | 85 | ~87 | ~39 | Change `set a 85.` in `setUpMix.tcl` |
| Base high conc. 3 | Neutral | Solvated | 74 | ~131 | ~50 | Change `set a 74.` in `setUpMix.tcl` |

## Prerequisites

- **NAMD 3.0** (CPU for equilibration, GPU for production)
- **VMD 1.9.3+** with `tempoUserVMD`, `tempoUtils`, `solvate`, `autoionize`, `psfgen` packages
- **R** with packages: `bio3d`, `doParallel`, `ggplot2`, `dplyr`, `tidyr`, `cowplot`, `Rcpp`, `sna`, `geometry`, `bigmemory`, `SOMMD`
- **CHARMM36 force field** (`par_all36_cgenff.prm`, `top_all36_cgenff.rtf`, `toppar_water_ions_jaf_cgenff.str`)

## Repository Structure

```
forcefield/         FFssFF-specific force field parameters and input PDBs
simulation/         Template VMD setup and NAMD config files
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

Run with NAMD (CPU multicore recommended):
```bash
namd3 +p <ncpus> npt_eq.inp
```

### 3. Production

Edit `simulation/npt01.inp`:
- Set `Dir` to match

Run with NAMD (GPU recommended):
```bash
namd3 +p <ncpus> +devices <gpu_id> npt01.inp
```

### 4. Vacuum Simulations

Use `nvt_eq.inp` and `nvt01.inp` instead of the NPT configs. No solvation step — use the lattice PSF/PDB directly.

## Analysis

### Contact Analysis

```bash
# For the 3 original systems (Acid, Base, Dry)
sbatch analysis/run_contact_analysis.slurm

# For the higher concentration sweep (high1, high2, high3)
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
| Cutoff | 12 Å |
| Thermostat | Langevin (γ = 1 ps⁻¹) |
| Barostat | Langevin piston (200 ps period, 50 ps decay) |
| Output frequency | 10,000 steps (20 ps) |

---

*README co-generated with LLM coding tools.*
