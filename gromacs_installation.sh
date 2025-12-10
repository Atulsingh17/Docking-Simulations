#!/bin/bash
# GROMACS 2024.4 CPU-only

# ... (same as before, just cmake line)
cmake .. -DGMX_GPU=OFF -DGMX_MPI=ON -DGMX_OPENMP=ON
make -j8
sudo make install
