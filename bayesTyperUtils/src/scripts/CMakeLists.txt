include_directories(${BayesTyperUtils_SOURCE_DIR}/vcf++/include/)
link_directories(${BayesTyperUtils_SOURCE_DIR}/vcf++/lib)

SET(EXECUTABLE_OUTPUT_PATH ${BayesTyperUtils_SOURCE_DIR}/bin/scripts)

add_executable(prepareSimulationFemaleSamples prepareSimulationFemaleSamples.cpp)
target_link_libraries(prepareSimulationFemaleSamples vcf++)

add_executable(getBiallelicVariantStats getBiallelicVariantStats.cpp)
target_link_libraries(getBiallelicVariantStats vcf++)

add_executable(getSimpleAlleleSummary getSimpleAlleleSummary.cpp)
target_link_libraries(getSimpleAlleleSummary vcf++)

add_executable(fixBayesTyperOutput fixBayesTyperOutput.cpp)
target_link_libraries(fixBayesTyperOutput vcf++)

add_executable(filterNotCalledAlternativeAlleles filterNotCalledAlternativeAlleles.cpp)
target_link_libraries(filterNotCalledAlternativeAlleles vcf++)

add_executable(splitMultiAllelicVariants splitMultiAllelicVariants.cpp)
target_link_libraries(splitMultiAllelicVariants vcf++)

add_executable(filterStructuralVariants filterStructuralVariants.cpp)
target_link_libraries(filterStructuralVariants vcf++)

add_executable(addUniqueIdToVariants addUniqueIdToVariants.cpp)
target_link_libraries(addUniqueIdToVariants vcf++)

add_executable(simplifyDoublePolymorphisms simplifyDoublePolymorphisms)
target_link_libraries(simplifyDoublePolymorphisms vcf++)

add_executable(getCallsetOverlapStats getCallsetOverlapStats)
target_link_libraries(getCallsetOverlapStats vcf++)

add_executable(evaluateCallset evaluateCallset)
target_link_libraries(evaluateCallset vcf++)

