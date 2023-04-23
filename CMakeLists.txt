cmake_minimum_required(VERSION 3.13)

project(planranger C CXX)
set(CMAKE_CXX_STANDARD 20)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

find_package(fmt 3.3 REQUIRED)
find_package(embree 3.0 REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Eigen3 3.3 REQUIRED)
find_package(ortools REQUIRED)

add_subdirectory(ext/)

add_executable(planranger
    src/main.cpp
    src/terrainmap.cpp
    src/terrainmesh.cpp
    src/config.cpp
    src/buildmaps.cpp
    src/path.cpp
    src/vantage.cpp
    src/plan.cpp
)

target_include_directories(planranger PRIVATE ext/ src/)
target_link_libraries(planranger fmt embree OpenMP::OpenMP_CXX miniz Eigen3::Eigen ortools::ortools)

if (UNIX AND NOT APPLE)
  target_link_libraries(planranger stdc++fs)
endif (UNIX AND NOT APPLE)

file(GLOB all_clang_files
  "src/*.h"
  "src/*.cpp")

add_custom_target(
        clang-format
        COMMAND /usr/bin/clang-format
        -style=file
        -i
        ${all_clang_files})

