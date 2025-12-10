#!/bin/bash
# AutoDock Vina + ADFR suite + Open Babel installation (Ubuntu/Debian/CentOS)

echo "Installing dependencies..."
sudo apt update && sudo apt install -y wget tar python3 python3-pip openbabel

# AutoDock Vina 1.2.5
wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64.bin
chmod +x vina_1.2.5_linux_x86_64.bin
sudo mv vina_1.2.5_linux_x86_64.bin /usr/local/bin/vina

# ADFRsuite (for receptor prep)
wget https://ccsb.scripps.edu/adfr/downloads/ADFRsuite_x86_64Linux_1.0.tar.gz
tar xzf ADFRsuite_x86_64Linux_1.0.tar.gz
cd ADFRsuite_x86_64Linux_1.0
sudo ./install.sh -d /usr/local -c 0

echo "Vina + ADFR installed! Test: vina --help"
