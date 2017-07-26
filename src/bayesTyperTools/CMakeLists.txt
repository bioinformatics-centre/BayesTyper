project(bayesTyperTools)

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
SET(EXTERNAL_DIR "${CMAKE_SOURCE_DIR}/external/")

include_directories(${CMAKE_SOURCE_DIR}/include/bayesTyperTools/ ${CMAKE_SOURCE_DIR}/include/vcf++ ${Boost_INCLUDE_DIRS} ${EXTERNAL_DIR})
link_directories(${CMAKE_SOURCE_DIR}/lib ${Boost_LIB_DIRS})

add_executable(${PROJECT_NAME} main.cpp convertAlleleId.cpp combine.cpp split.cpp merge.cpp annotate.cpp writeIndels.cpp addAttributes.cpp filter.cpp getSummary.cpp ${EXTERNAL_DIR}/edlib/edlib.cpp)

target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES} vcf++)
