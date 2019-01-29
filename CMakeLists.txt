cmake_minimum_required(VERSION 2.8.0)


execute_process(
  COMMAND git rev-parse HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_LAST_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)


SET(KMER_SIZE 55 CACHE STRING "BayesTyper kmer size")

if(${KMER_SIZE} GREATER 255)

	message(FATAL_ERROR "\nKmer size needs to be <= 255")

else(${KMER_SIZE} GREATER 255)

	message(STATUS "Kmer size: " ${KMER_SIZE})

endif(${KMER_SIZE} GREATER 255)


SET(BUILD_STATIC 0 CACHE BOOL "Build static")

if(${BUILD_STATIC} EQUAL 0)

	SET(CMAKE_CXX_FLAGS "--std=c++11 -lpthread -g -O3 -DBT_KMER_SIZE='${KMER_SIZE}' -DBT_VERSION='\"v1.4.1 ${GIT_LAST_COMMIT_HASH}\"'")

else(${BUILD_STATIC} EQUAL 0)

	message(STATUS "Building BayesTyper static")
	
	SET(CMAKE_CXX_FLAGS "--std=c++11 -pthread -O3 -static -static-libgcc -static-libstdc++ -DBT_KMER_SIZE='${KMER_SIZE}' -DBT_VERSION='\"v1.4.1\"'")
	SET(CMAKE_EXE_LINKER_FLAGS "-Wl,--whole-archive -lpthread -Wl,--no-whole-archive")

	SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
	SET(BUILD_SHARED_LIBS OFF)
	SET(Boost_USE_STATIC_LIBS ON)

endif(${BUILD_STATIC} EQUAL 0)


FIND_PACKAGE(Boost COMPONENTS iostreams program_options system filesystem serialization REQUIRED)
message(STATUS ${Boost_LIBRARIES})

add_subdirectory(${CMAKE_SOURCE_DIR}/external/kmc_api)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/vcf++)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/kmerBloom)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/bayesTyper)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/bayesTyperTools)


SET(BUILD_SCRIPTS 0 CACHE BOOL "Build BayesTyperTools scripts")

if(${BUILD_SCRIPTS} EQUAL 1)

	add_subdirectory(${CMAKE_SOURCE_DIR}/src/bayesTyperTools/scripts)

endif(${BUILD_SCRIPTS} EQUAL 1)

