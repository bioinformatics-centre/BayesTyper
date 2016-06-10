cmake_minimum_required (VERSION 2.8.0)

set(CMAKE_CXX_COMPILER "c++")

project(BayesTyperUtils)

set(CMAKE_CXX_FLAGS "--std=c++11 -g -Wall -O3 -DBTU_VERSION='\"v0.9\"' -lpthread")

add_subdirectory(vcf++)
add_subdirectory(src)
add_subdirectory(src/scripts)
