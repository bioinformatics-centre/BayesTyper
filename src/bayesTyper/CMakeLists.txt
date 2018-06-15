
project(bayesTyper)

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
	
include_directories(${CMAKE_SOURCE_DIR}/include/bayesTyper ${CMAKE_SOURCE_DIR}/include/vcf++ ${CMAKE_SOURCE_DIR}/include/kmerBloom ${CMAKE_SOURCE_DIR}/external ${Boost_INCLUDE_DIRS})
link_directories(${CMAKE_SOURCE_DIR}/lib)

add_executable(${PROJECT_NAME} main.cpp VariantClusterGenotyper.cpp VariantClusterGroup.cpp InferenceEngine.cpp CountDistribution.cpp FrequencyDistribution.cpp SparsityEstimator.cpp Sample.cpp DiscreteSampler.cpp Combinator.cpp GenotypeWriter.cpp CountAllocation.cpp HaplotypeFrequencyDistribution.cpp KmerCounts.cpp KmerCounter.cpp KmerHash.cpp VariantClusterGraph.cpp VariantFileParser.cpp VariantClusterHaplotypes.cpp ChromosomePloidy.cpp KmerStats.cpp NegativeBinomialDistribution.cpp Regions.cpp Chromosomes.cpp VariantClusterGraphPath.cpp Filters.cpp)

target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES} vcf++ kmc kmerBloom z)
