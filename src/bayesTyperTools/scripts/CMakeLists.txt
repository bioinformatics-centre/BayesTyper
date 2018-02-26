
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin/scripts)

include_directories(${CMAKE_SOURCE_DIR}/include/bayesTyperTools ${CMAKE_SOURCE_DIR}/include/vcf++ ${CMAKE_SOURCE_DIR}/external/ ${Boost_INCLUDE_DIRS})
link_directories(${CMAKE_SOURCE_DIR}/lib ${Boost_LIB_DIRS})

add_executable(generateDiplotypes generateDiplotypes.cpp)
target_link_libraries(generateDiplotypes vcf++)

add_executable(filterStructuralVariants filterStructuralVariants.cpp)
target_link_libraries(filterStructuralVariants vcf++)

add_executable(getKmerStats getKmerStats.cpp)
target_link_libraries(getKmerStats vcf++ kmc)

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

add_executable(transferAlleleCallsetOrigin transferAlleleCallsetOrigin.cpp)
target_link_libraries(transferAlleleCallsetOrigin vcf++)

add_executable(writeIndels writeIndels.cpp)
target_link_libraries(writeIndels vcf++)


