# Create AsciiReader library.
FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)
#list(REMOVE_ITEM lib_srcs "${CMAKE_CURRENT_SOURCE_DIR}/DataTest.cpp")
add_library (AsciiReader SHARED ${lib_srcs} ${lib_headers})
target_link_libraries( AsciiReader DataReader )

target_include_directories (AsciiReader PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

#install (TARGETS AsciiReader DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)
