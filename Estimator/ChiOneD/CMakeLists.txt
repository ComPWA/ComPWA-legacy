# Create ChiOneD library.
FILE(GLOB lib_srcs *.cpp)
add_library (ChiOneD SHARED ${lib_srcs})
target_link_libraries( ChiOneD )

target_include_directories (ChiOneD PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})