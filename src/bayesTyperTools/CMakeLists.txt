project(bayesTyperTools)

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

include_directories(${CMAKE_SOURCE_DIR}/include/bayesTyperTools ${CMAKE_SOURCE_DIR}/include/vcf++ ${CMAKE_SOURCE_DIR}/include/kmerBloom ${CMAKE_SOURCE_DIR}/external ${Boost_INCLUDE_DIRS})
link_directories(${CMAKE_SOURCE_DIR}/lib ${Boost_LIB_DIRS})

add_executable(${PROJECT_NAME} main.cpp MakeBloom.cpp ConvertAllele.cpp Combine.cpp Annotate.cpp AddAttributes.cpp Filter.cpp ${CMAKE_SOURCE_DIR}/external/edlib/edlib.cpp)

target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES} vcf++ kmc kmerBloom)
