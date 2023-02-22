cmake_minimum_required(VERSION 3.13)

add_subdirectory(tomlplusplus)

add_library(miniz miniz.c miniz.h)
target_include_directories(miniz PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
