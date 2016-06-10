#SET(EXECUTABLE_OUTPUT_PATH ${VCF++_SOURCE_DIR}/bin)
include_directories(${VCF++_SOURCE_DIR}/include/vcf++)
link_directories(${VCF++_SOURCE_DIR}/lib)

add_executable(runUnitTests catch.cpp StatsTest.cpp CompareOperatorsTest.cpp)
target_link_libraries(runUnitTests vcf++)
