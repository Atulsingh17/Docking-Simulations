# Computational Drug Discovery Toolkit  
**Molecular Docking · MD Simulations · Binding Free Energy · Trajectory Analysis**

PhD-level in-silico drug discovery pipelines  
Atul Kumar Singh | Computational & Molecular Biology  
GitHub: https://github.com/Atulsingh17

## Quick Overview
| Repository / Tool              | Purpose                                      | Key Software          |
|--------------------------------|----------------------------------------------|------------------------|
| Molecular Docking              | High-throughput virtual screening            | AutoDock Vina + Gnina |
| GROMACS Simulations            | All-atom & GPU-accelerated MD                | GROMACS 2024.4 (GPU)  |
| MM/GBSA & PBSA                 | Accurate binding free energy calculation     | gmx_MMPBSA (Amber)    |
| Complete Trajectory Analysis   | RMSD/RMSF/SASA/Rg/PCA/Free Energy Landscape  | MDAnalysis + Python   |

## One-Command Installations

```bash
# AutoDock Vina + ADFR suite
bash install_autodock_vina.sh

# GROMACS 2024.4 with GPU support
bash install_gromacs_gpu.sh

# GROMACS CPU-only
bash install_gromacs_cpu.sh
