# Create ComPWA Core library.
FILE(GLOB lib_srcs *.cpp)
#FILE(GLOB lib_ignore *Test.cpp)
list(REMOVE_ITEM lib_srcs "${CMAKE_CURRENT_SOURCE_DIR}/ParameterListTest.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/ParameterTest.cpp")
add_library (Core SHARED ${lib_srcs})
target_link_libraries( Core ${Boost_LIBRARIES} )

target_include_directories (Core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})