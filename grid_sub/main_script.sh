#!/bin/bash

pwd
echo "----------"
ls -altrh /
echo "----------"
ls `pwd`

path=`pwd`

#cat /cvmfs/sft.cern.ch/lcg/views/LCG_98/x86_64-slc6-gcc8-opt/setup.sh
cat /cvmfs/sft.cern.ch/lcg/views/LCG_96/x86_64-slc6-gcc8-opt/setup.sh

source /cvmfs/sft.cern.ch/lcg/views/LCG_96/x86_64-slc6-gcc8-opt/setup.sh
#source /cvmfs/sft.cern.ch/lcg/views/LCG_98/x86_64-slc6-gcc8-opt/setup.sh


#ls /cvmfs
ls /cvmfs/sft.cern.ch/lcg/views/LCG_96/x86_64-slc6-gcc8-opt
echo "+++++++"
ls /cvmfs/sft.cern.ch/lcg/views/LCG_96/
echo "-------"
ls /cvmfs/sft.cern.ch/lcg/views/LCG_98/
echo "======="
ls /cvmfs/sft.cern.ch/lcg/views

#git clone -b feature/nova_submodules git@github.com:P-theta/P-theta.git P-theta

mkdir P-theta_install

cd P-theta
#git submodule update --init --recursive
#cat $path/P-theta/Minimal/library/applications/PTBestFit.h
[ -d build ] || rm -rf build
mkdir build
cd build

#-DCMAKE_BUILD_TYPE=DEBUG
cmake -DCMAKE_BUILD_TYPE=DEBUG -D CMAKE_CXX_STANDARD=17 -D PYBIND11_INSTALL=FALSE -D USE_NOVATOOLS=1 -D Scarab_BUILD_PYTHON=FALSE -D CMAKE_INSTALL_PREFIX=$path/P-theta_install ../Minimal
cmake -DCMAKE_BUILD_TYPE=DEBUG -D CMAKE_CXX_STANDARD=17 -D PYBIND11_INSTALL=FALSE -D USE_NOVATOOLS=1 -D Scarab_BUILD_PYTHON=FALSE -D CMAKE_INSTALL_PREFIX=$path/P-theta_install ../Minimal
#cmake -D PYBIND11_INSTALL=FALSE -D USE_NOVATOOLS=1 -D Scarab_BUILD_PYTHON=FALSE -D CMAKE_INSTALL_PREFIX=../install  ../Minimal
make -j3
make -j3 install

echo "Is the installation successful?"
ls $path/P-theta_install
ls $path/P-theta_install/bin
ls $path/P-theta_install/include
ls $path/P-theta_install/lib

echo "**cat CMakeCache.txt"
cat CMakeCache.txt
echo "**cat library/externals/scarab/CMakeLists.txt"
cat ../Minimal/library/externals/scarab/CMakeLists.txt
echo "**cat library/externals/mcmc/CMakeLists.txt"
cat ../Minimal/library/externals/mcmc/CMakeLists.txt
echo "**cat library/externals/ModProb3++/CMakeLists.txt"
cat ../Minimal/library/externals/ModProb3++/CMakeLists.txt

cd ../../

tar -czvf P-theta_install.gzip P-theta_install
