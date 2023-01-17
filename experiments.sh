#!/usr/bin/env sh

rm -rf results/

# Slope Experiment
mkdir -p results/slope05
./build/planranger -p=0.25 -va=55 -nv=15 -rs=5 meshes/lmp.ply results/slope05/

mkdir -p results/slope10
./build/planranger -p=0.25 -va=55 -nv=15 -rs=10 meshes/lmp.ply results/slope10/

mkdir -p results/slope15
./build/planranger -p=0.25 -va=55 -nv=15 -rs=15 meshes/lmp.ply results/slope15/

mkdir -p results/slope20
./build/planranger -p=0.25 -va=55 -nv=15 -rs=20 meshes/lmp.ply results/slope20/
