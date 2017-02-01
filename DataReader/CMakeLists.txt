# Create DataReader library.
FILE(GLOB lib_srcs *.cpp)
add_library (DataReader ${lib_srcs})

add_subdirectory (RootReader)
add_subdirectory (JakeReader)
#add_subdirectory (AsciiReader)

target_link_libraries( DataReader ${Boost_LIBRARIES} RootReader JakeReader)

target_include_directories (DataReader PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})