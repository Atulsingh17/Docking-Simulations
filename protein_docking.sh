#!/bin/bash
# Full docking pipeline for protein (PDB: e.g., 6E2F example)

# Step 1: Prepare receptor
prepare_receptor -r protein.pdb -o protein.pdbqt

# Step 2: Prepare ligand (SMILES or SDF)
obabel ligand.sdf -O ligand.pdbqt --gen3d

# Step 3: Define grid box (protein binding pocket)
cat > config.txt << EOF
center_x = 00.0
center_y = 00.0
center_z = 00.0
size_x = 00
size_y = 00
size_z = 00
exhaustiveness = 00
EOF

# Step 4: Run docking
vina --receptor protein.pdbqt --ligand ligand.pdbqt --config config.txt --out docking_output.pdbqt --log docking_log.txt

echo "Docking complete! Best pose:"
grep "1 " docking_log.txt
