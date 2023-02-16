# PitMissionPlanner

## Build

Clone the code!
```
git clone --recursive git@github.com:PlanetaryRobotics/PitMissionPlanner.git
```

### Ubuntu 22.04
Install build tools and dependencies.
```
sudo apt install libfmt-dev libeigen3-dev libembree3-dev cmake g++
```

Create a build directory and cd into it. Then, run cmake and make to build the code.
```
mkdir build
cd build
cmake ..
make
```

### M1 Mac
Install build tools and dependencies.
```
brew install fmt embree eigen3 llvm cmake
```

Create a build directory and cd into it. Then, run cmake and make to build the code.
```
mkdir build
cd build
cmake ..
make
```

### Windows
Enable windows subsystem for linux and install the Ubuntu distro. Then follow the Ubuntu build instructions!

## Viewing the Results
The .exr maps produced by the planner can be viewed using the [tev image viewer](https://github.com/Tom94/tev).

The .ply meshes and .xyz point clouds can be viewed using [CloudCompare](https://github.com/CloudCompare/CloudCompare).
