SET(LIBRARY_OUTPUT_PATH ${VCF++_SOURCE_DIR}/lib)

include_directories(${VCF++_SOURCE_DIR}/include/vcf++)
add_library(vcf++ VcfFile.cpp Variant.cpp Allele.cpp VcfMetaData.cpp AttributeSet.cpp Sample.cpp Utils.cpp Attribute.cpp JoiningString.cpp Stats.cpp AttributeFilter.cpp SampleAlleleAttributeFilter.cpp FastaReader.cpp FastaRecord.cpp Auxiliaries.cpp Trio.cpp Contig.cpp Regions.cpp)
