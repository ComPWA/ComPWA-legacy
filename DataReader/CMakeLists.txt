# Create DataReader library.
FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)
add_library (ComPWA_DataReader SHARED ${lib_srcs} ${lib_headers})

add_subdirectory (RootReader)
add_subdirectory (JakeReader)
add_subdirectory (AsciiReader)

target_link_libraries( ComPWA_DataReader ComPWA_Core ${Boost_LIBRARIES} )

target_include_directories (ComPWA_DataReader PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

#install (TARGETS ComPWA_DataReader DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)
