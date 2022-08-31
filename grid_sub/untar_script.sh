#!/bin/bash
echo 'Untar step:'
ls -atr
python -c 'import tarfile; tar = tarfile.open("Packages.gzip"); tar.extractall(); tar.close()'
python -c 'import tarfile; tar = tarfile.open("inputs_Spring2020_0325.gzip"); tar.extractall(); tar.close()'
python -c 'import tarfile; tar = tarfile.open("nova-sl7-novat2k_v4_fixcosmicsrock.tar.bz2"); tar.extractall(); tar.close()'
echo 'Untar step done'
