#!/usr/bin/env sh

rm -rf results/

# Slope Experiment
mkdir -p results/slope05
./build/planranger -p=1.0 -va=55 -nv=18 -lonslope=5  -latslope=5  -turnslope=8 meshes/lmp_crop.ply results/slope05/

mkdir -p results/slope10
./build/planranger -p=1.0 -va=55 -nv=18 -lonslope=10 -latslope=8  -turnslope=8 meshes/lmp_crop.ply results/slope10/

mkdir -p results/slope15
./build/planranger -p=1.0 -va=55 -nv=18 -lonslope=15 -latslope=10 -turnslope=8 meshes/lmp_crop.ply results/slope15/

mkdir -p results/slope20
./build/planranger -p=1.0 -va=55 -nv=18 -lonslope=20 -latslope=10 -turnslope=8 meshes/lmp_crop.ply results/slope20/
