project(vcf++)

SET(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

include_directories(${CMAKE_SOURCE_DIR}/include/vcf++ ${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIB_DIRS})

add_library(${PROJECT_NAME} VcfFile.cpp Variant.cpp Allele.cpp VcfMetaData.cpp AttributeSet.cpp Sample.cpp Utils.cpp Attribute.cpp JoiningString.cpp Stats.cpp AttributeFilter.cpp SampleAlleleAttributeFilter.cpp FastaReader.cpp FastaRecord.cpp Auxiliaries.cpp Trio.cpp Contig.cpp Regions.cpp)
