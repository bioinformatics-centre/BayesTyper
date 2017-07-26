
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin/scripts)
SET(EXTERNAL_DIR "${CMAKE_SOURCE_DIR}/external/")

include_directories(${CMAKE_SOURCE_DIR}/include/bayesTyperTools/ ${CMAKE_SOURCE_DIR}/include/vcf++ ${Boost_INCLUDE_DIRS} ${EXTERNAL_DIR})
link_directories(${CMAKE_SOURCE_DIR}/lib ${Boost_LIB_DIRS})

add_executable(generateDiplotypes generateDiplotypes.cpp)
target_link_libraries(generateDiplotypes vcf++)

add_executable(filterStructuralVariants filterStructuralVariants.cpp)
target_link_libraries(filterStructuralVariants vcf++)

add_executable(getKmerStats getKmerStats.cpp ${EXTERNAL_DIR}/kmc_api/kmc_file.cpp ${EXTERNAL_DIR}/kmc_api/kmer_api.cpp ${EXTERNAL_DIR}/kmc_api/mmer.cpp)
target_link_libraries(getKmerStats vcf++)

add_executable(transferAlleleFreqFields transferAlleleFreqFields.cpp)
target_link_libraries(transferAlleleFreqFields vcf++)

add_executable(collapseSummaryTable collapseSummaryTable.cpp)
target_link_libraries(collapseSummaryTable vcf++)

add_executable(assessHaplotypeTransmissionSupport assessHaplotypeTransmissionSupport.cpp)
target_link_libraries(assessHaplotypeTransmissionSupport vcf++)

add_executable(getGenomicIntervals getGenomicIntervals.cpp)
target_link_libraries(getGenomicIntervals vcf++)

add_executable(convertNestedGenotypes convertNestedGenotypes.cpp)
target_link_libraries(convertNestedGenotypes vcf++)

add_executable(addMaxGenotypePosterior addMaxGenotypePosterior.cpp)
target_link_libraries(addMaxGenotypePosterior vcf++)

add_executable(addEditDistanceAndCondordance addEditDistanceAndCondordance.cpp)
target_link_libraries(addEditDistanceAndCondordance vcf++)

add_executable(convertSeqToAlleleId convertSeqToAlleleId.cpp)
target_link_libraries(convertSeqToAlleleId vcf++)

add_executable(addGenotypeQuality addGenotypeQuality.cpp)
target_link_libraries(addGenotypeQuality vcf++)



