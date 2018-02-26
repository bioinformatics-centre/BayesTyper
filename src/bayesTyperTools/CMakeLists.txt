project(bayesTyperTools)

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

include_directories(${CMAKE_SOURCE_DIR}/include/bayesTyperTools ${CMAKE_SOURCE_DIR}/include/vcf++ ${CMAKE_SOURCE_DIR}/include/kmerBloom ${CMAKE_SOURCE_DIR}/external ${CMAKE_SOURCE_DIR}/external/libbf ${Boost_INCLUDE_DIRS})
link_directories(${CMAKE_SOURCE_DIR}/lib ${Boost_LIB_DIRS})

add_executable(${PROJECT_NAME} main.cpp bloomMaker.cpp convertAlleleId.cpp combine.cpp split.cpp merge.cpp annotate.cpp addAttributes.cpp filter.cpp getSummary.cpp ${CMAKE_SOURCE_DIR}/external/edlib/edlib.cpp)

target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES} vcf++ kmc kmerBloom libbf)
