
project(bayesTyper)

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSMALLMERSIZE=18")

SET(EXTERNAL_DIR "${CMAKE_SOURCE_DIR}/external/")

include_directories(${CMAKE_SOURCE_DIR}/include/bayesTyper/ ${CMAKE_SOURCE_DIR}/include/vcf++ ${Boost_INCLUDE_DIRS} ${EXTERNAL_DIR})
link_directories(${CMAKE_SOURCE_DIR}/lib ${Boost_LIB_DIRS})

add_executable(${PROJECT_NAME} main.cpp VariantClusterGenotyper.cpp VariantClusterGroup.cpp InferenceEngine.cpp CountDistribution.cpp FrequencyDistribution.cpp SparsityEstimator.cpp DiscreteSampler.cpp Combinator.cpp GenotypeWriter.cpp CountAllocation.cpp VariantCluster.cpp HaplotypeFrequencyDistribution.cpp KmerCounts.cpp KmerCounter.cpp KmerHash.cpp VariantClusterGraph.cpp VariantFileParser.cpp KmerFactory.cpp VariantClusterHaplotypes.cpp ChromosomePloidy.cpp KmerStats.cpp NegativeBinomialDistribution.cpp Regions.cpp VariantClusterGraphPath.cpp ${EXTERNAL_DIR}/kmc_api/kmc_file.cpp ${EXTERNAL_DIR}/kmc_api/kmer_api.cpp ${EXTERNAL_DIR}/kmc_api/mmer.cpp)

target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES} vcf++)
