# Create DataReader library.
FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)
add_library (DataReader SHARED ${lib_srcs} ${lib_headers})

add_subdirectory (RootReader)
add_subdirectory (JakeReader)
add_subdirectory (AsciiReader)

target_link_libraries( DataReader Core ${Boost_LIBRARIES} )

target_include_directories (DataReader PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

#install (TARGETS DataReader DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)
