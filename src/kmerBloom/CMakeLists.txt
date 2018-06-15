project(kmerBloom)

SET(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

include_directories(${CMAKE_SOURCE_DIR}/include/kmerBloom ${CMAKE_SOURCE_DIR}/external/)
link_directories(${CMAKE_SOURCE_DIR}/lib)

add_library(${PROJECT_NAME} KmerBloom.cpp)
