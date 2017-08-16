cmake_minimum_required(VERSION 2.8.0)

execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
  COMMAND git rev-parse HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_LAST_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

SET(CMAKE_CXX_FLAGS "--std=c++11 -g -Wall -O3 -DBT_VERSION='\"v1.1 ${GIT_LAST_COMMIT_HASH}\"' -lpthread")

FIND_PACKAGE(Boost COMPONENTS program_options system filesystem iostreams REQUIRED)
message(STATUS ${Boost_LIBRARIES})
add_subdirectory(${CMAKE_SOURCE_DIR}/src/vcf++)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/bayesTyper)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/bayesTyperTools)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/bayesTyperTools/scripts)
