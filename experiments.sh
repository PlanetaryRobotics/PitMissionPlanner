#!/usr/bin/env sh

export OMP_NUM_THREADS=4

rm -rf results/

# Slope Experiments
mkdir -p results/slope05
./build/planranger -lon=5 --OutputDir=results/slope05/ ./missions/interior.toml
./scripts/makevideo.sh results/slope05/anim/

mkdir -p results/slope10
./build/planranger -lon=10 --OutputDir=results/slope10/ ./missions/interior.toml
./scripts/makevideo.sh results/slope10/anim/

mkdir -p results/slope15
./build/planranger -lon=15 --OutputDir=results/slope15/ ./missions/interior.toml
./scripts/makevideo.sh results/slope15/anim/

mkdir -p results/slope20
./build/planranger -lon=20 --OutputDir=results/slope20/ ./missions/interior.toml
./scripts/makevideo.sh results/slope20/anim/

# Science Priority Experiments
mkdir -p results/science-interior
./build/planranger --OutputDir=results/science-interior/ ./missions/interior.toml
./scripts/makevideo.sh results/science-interior/anim/

mkdir -p results/science-notch
./build/planranger --OutputDir=results/science-notch/ ./missions/notch.toml
./scripts/makevideo.sh results/science-notch/anim/

mkdir -p results/science-ramp
./build/planranger --OutputDir=results/science-ramp/ ./missions/ramp.toml
./scripts/makevideo.sh results/science-ramp/anim/

mkdir -p results/science-cave
./build/planranger --OutputDir=results/science-cave/ ./missions/cave.toml
./scripts/makevideo.sh results/science-cave/anim/
