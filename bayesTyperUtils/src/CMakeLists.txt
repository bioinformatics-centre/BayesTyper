find_package(Boost COMPONENTS program_options system filesystem)

include_directories(${BayesTyperUtils_SOURCE_DIR}/vcf++/include/ ${Boost_INCLUDE_DIRS})
link_directories(${BayesTyperUtils_SOURCE_DIR}/vcf++/lib ${Boost_LIBRARY_DIRS})

SET(EXECUTABLE_OUTPUT_PATH ${BayesTyperUtils_SOURCE_DIR}/bin)
add_executable(BayesTyperUtils main.cpp convertAlleleIdToSequence.cpp clean.cpp remove.cpp combine.cpp split.cpp merge.cpp annotate.cpp calcFilterStats.cpp filter.cpp markDeNovo.cpp writeIndels.cpp addRepeatAnnotation.cpp generateRelease.cpp getSummary.cpp calcTrioConcordance.cpp selectValidationVariants.cpp)

target_link_libraries(BayesTyperUtils ${Boost_LIBRARIES} vcf++)
