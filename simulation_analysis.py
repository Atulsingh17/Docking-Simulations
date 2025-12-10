# Comprehensive MD simulation analysis: RMSD, RMSF, SASA, Rg, PCA, FEL
# Requirements: pip install MDAnalysis numpy matplotlib scikit-learn
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align, pca, gyrate, sasa
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.stats import binned_statistic_2d

# Load universe (replace with your files)
u = mda.Universe("md_100ns.tpr", "md_100ns.xtc")
protein = u.select_atoms("protein")
backbone = u.select_atoms("backbone")

# 1. RMSD
align.AlignTraj(u, u, select="backbone", in_memory=True).run()
rmsd = rms.RMSD(u, select="backbone").run()
plt.figure()
plt.plot(rmsd.rmsd[:, 1], rmsd.rmsd[:, 2])
plt.xlabel("Time (ps)"); plt.ylabel("RMSD (Å)")
plt.title("RMSD of Protein Backbone")
plt.savefig("rmsd.png")

# 2. RMSF
rmsf = rms.RMSF(backbone).run()
plt.figure()
plt.plot(rmsf.rmsf)
plt.xlabel("Residue"); plt.ylabel("RMSF (Å)")
plt.title("RMSF of Protein Backbone")
plt.savefig("rmsf.png")

# 3. SASA (Solvent Accessible Surface Area)
sasa_analysis = sasa.SASA(u, probe_radius=1.4, n_slices=100)
sasa_analysis.run()
plt.figure()
plt.plot(sasa_analysis.times, sasa_analysis.total_area)
plt.xlabel("Time (ps)"); plt.ylabel("SASA (Å²)")
plt.title("Total SASA of Protein")
plt.savefig("sasa.png")

# 4. Radius of Gyration (Rg)
rg = gyrate.Gyrate(protein)
rg.run()
plt.figure()
plt.plot(rg.times, rg.rg)
plt.xlabel("Time (ps)"); plt.ylabel("Rg (Å)")
plt.title("Radius of Gyration of Protein")
plt.savefig("rg.png")

# 5. PCA (Principal Component Analysis)
pc = pca.PCA(u, select="backbone")
pc.run()
variance = pc.cumulated_variance[:10]  # Top 10 modes
plt.figure()
plt.bar(range(1, 11), variance)
plt.xlabel("Principal Component"); plt.ylabel("Cumulative Variance")
plt.title("PCA Cumulative Variance")
plt.savefig("pca_variance.png")

# Project trajectory on PC1 and PC2
atomgroup = u.select_atoms("backbone")
transformed = pc.transform(atomgroup, n_components=2)
pc1, pc2 = transformed[:, 0], transformed[:, 1]

# 6. Free Energy Landscape (FEL) from PC1/PC2
delta_g, bins_pc1, bins_pc2 = binned_statistic_2d(pc1, pc2, None, statistic='count', bins=50)
delta_g = -np.log(delta_g / np.max(delta_g))  # kT = 1 assumption
delta_g[np.isinf(delta_g)] = np.nan

plt.figure()
extent = [bins_pc1.min(), bins_pc1.max(), bins_pc2.min(), bins_pc2.max()]
plt.imshow(delta_g.T, origin='lower', extent=extent, cmap='viridis')
plt.colorbar(label="ΔG (kT)")
plt.xlabel("PC1"); plt.ylabel("PC2")
plt.title("Free Energy Landscape (FEL)")
plt.savefig("fel.png")

print("All analyses complete! Check .png files for plots.")
