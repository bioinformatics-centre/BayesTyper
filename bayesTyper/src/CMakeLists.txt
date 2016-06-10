cmake_minimum_required(VERSION 2.6)

project(bayesTyper)

SET(KMC_API_DIR "../kmc_api/")
SET(CMAKE_CXX_COMPILER "g++")

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSMALLMERSIZE=18")

FIND_PACKAGE(Boost COMPONENTS program_options iostreams REQUIRED)
message(STATUS ${Boost_LIBRARIES})

include_directories(. ${Boost_INCLUDE_DIRS} ../ProducerConsumerQueue/src/ ${KMC_API_DIR})
link_directories(${Boost_LIB_DIRS})

add_executable(${PROJECT_NAME} main.cpp VariantClusterGenotyper.cpp VariantClusterGroup.cpp InferenceEngine.cpp CountDistribution.cpp FrequencyDistribution.cpp SparsityEstimator.cpp DiscreteSampler.cpp Combinator.cpp GenotypeWriter.cpp CountAllocation.cpp VariantCluster.cpp HaplotypeFrequencyDistribution.cpp KmerCounts.cpp KmerCounter.cpp KmerHash.cpp VariantClusterGraph.cpp VariantFileParser.cpp KmerFactory.cpp VariantClusterHaplotypes.cpp ChromosomePloidy.cpp VariantKmerStats.cpp VariantClusterGraphVertex.cpp KmerCoverage.cpp NegativeBinomialDistribution.cpp Regions.cpp PartitionGraph.cpp ${KMC_API_DIR}kmc_file.cpp ${KMC_API_DIR}kmer_api.cpp ${KMC_API_DIR}mmer.cpp)

#target_link_libraries(${PROJECT_NAME} ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_IOSTREAMS_LIBRARY})
target_link_libraries(${PROJECT_NAME} boost_program_options boost_iostreams pthread)
