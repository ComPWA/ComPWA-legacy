# Create ComPWA Core library.
FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)
#FILE(GLOB lib_ignore *Test.cpp)
list(REMOVE_ITEM lib_srcs "${CMAKE_CURRENT_SOURCE_DIR}/ParameterListTest.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/ParameterTest.cpp")

add_library ( Core SHARED ${lib_srcs} ${lib_headers} )

target_link_libraries( Core ${Boost_LIBRARIES} )

target_include_directories ( Core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

#install (TARGETS Core DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)
