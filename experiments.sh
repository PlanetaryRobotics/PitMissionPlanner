#!/usr/bin/env sh

export OMP_NUM_THREADS=12

rm -rf results/

# Landing Site Experiments
mkdir -p results/site-good
./build/planranger -x=515 -y=395 -lh=2 --OutputDir=results/site-good/ ./missions/interior.toml | tee ./results/site-good/log.txt
./scripts/makevideo.sh results/site-good/anim/

mkdir -p results/site-bad
./build/planranger -x=480 -y=360 -lh=2 --OutputDir=results/site-bad/ ./missions/interior.toml | tee ./results/site-bad/log.txt
./scripts/makevideo.sh results/site-bad/anim/

# Slope Experiments
mkdir -p results/slope05
./build/planranger -lon=5 --OutputDir=results/slope05/ ./missions/interior.toml | tee ./results/slope05/log.txt
./scripts/makevideo.sh results/slope05/anim/

mkdir -p results/slope10
./build/planranger -lon=10 --OutputDir=results/slope10/ ./missions/interior.toml | tee ./results/slope10/log.txt
./scripts/makevideo.sh results/slope10/anim/

mkdir -p results/slope15
./build/planranger -lon=15 --OutputDir=results/slope15/ ./missions/interior.toml | tee ./results/slope15/log.txt
./scripts/makevideo.sh results/slope15/anim/

mkdir -p results/slope20
./build/planranger -lon=20 --OutputDir=results/slope20/ ./missions/interior.toml | tee ./results/slope20/log.txt
./scripts/makevideo.sh results/slope20/anim/

# Science Priority Experiments
mkdir -p results/science-interior
./build/planranger --OutputDir=results/science-interior/ ./missions/interior.toml | tee ./results/science-interior/log.txt
./scripts/makevideo.sh results/science-interior/anim/

mkdir -p results/science-notch
./build/planranger --OutputDir=results/science-notch/ ./missions/notch.toml | tee ./results/science-notch/log.txt
./scripts/makevideo.sh results/science-notch/anim/

mkdir -p results/science-ramp
./build/planranger --OutputDir=results/science-ramp/ ./missions/ramp.toml | tee ./results/science-ramp/log.txt
./scripts/makevideo.sh results/science-ramp/anim/

mkdir -p results/science-cave
./build/planranger --OutputDir=results/science-cave/ ./missions/cave.toml | tee ./results/science-cave/log.txt
./scripts/makevideo.sh results/science-cave/anim/
