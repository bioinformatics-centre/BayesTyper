# -- Project Setup ------------------------------------------------------------

project(libbf)

SET(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

include_directories(${CMAKE_SOURCE_DIR}/external/libbf)

set(libbf_sources
  src/bitvector.cpp
  src/counter_vector.cpp
  src/hash.cpp
  src/bloom_filter/a2.cpp
  src/bloom_filter/basic.cpp
  src/bloom_filter/bitwise.cpp
  src/bloom_filter/counting.cpp
  src/bloom_filter/stable.cpp
)

add_library(libbf ${libbf_sources})
