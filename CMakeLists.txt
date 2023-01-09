cmake_minimum_required(VERSION 3.13)

project(planranger CXX)
set(CMAKE_CXX_STANDARD 20)

find_package(fmt 3.3 REQUIRED)
find_package(embree 3.0 REQUIRED)
find_package(OpenMP REQUIRED)

add_executable(planranger src/main.cpp src/terrainmesh.cpp)
target_link_libraries(planranger fmt stdc++fs embree OpenMP::OpenMP_CXX)
