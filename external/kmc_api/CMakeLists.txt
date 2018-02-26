project(kmc)

SET(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

include_directories(${CMAKE_SOURCE_DIR}/external/kmc_api)

add_library(${PROJECT_NAME} kmc_file.cpp kmer_api.cpp mmer.cpp)
