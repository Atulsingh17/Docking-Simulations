#!/bin/bash
# GROMACS 2024.4 with CUDA 12 (NVIDIA GPU)

wget https://ftp.gromacs.org/gromacs/gromacs-2024.4.tar.gz
tar xzf gromacs-2024.4.tar.gz
cd gromacs-2024.4
mkdir build-gpu && cd build-gpu

cmake .. -DGMX_BUILD_OWN_FFT=ON -DREGRESSIONTEST_DOWNLOAD=ON -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-12
make -j8
sudo make install

source /usr/local/gromacs/bin/GMXRC
echo "GROMACS GPU ready! Test: gmx --version"
