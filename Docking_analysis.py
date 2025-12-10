# Docking result analysis + best pose extraction
import os
from openbabel import openbabel as ob

def extract_best_pose(pdbqt_file, output_pdb="best_pose.pdb"):
    poses = []
    with open(pdbqt_file) as f:
        lines = f.readlines()
    model = []
    for line in lines:
        if line.startswith("MODEL"):
            if model: poses.append("\n".join(model))
            model = [line]
        elif line.startswith("ENDMDL"):
            model.append(line)
    if model: poses.append("\n".join(model))
    
    # Save best pose (MODEL 1)
    with open(output_pdb, "w") as f:
        f.write(poses[0].replace("MODEL        1", ""))
    print(f"Best pose saved as {output_pdb}")

extract_best_pose("docking_output.pdbqt")
