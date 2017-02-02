# Create MinLogLH library.
FILE(GLOB lib_srcs *.cpp)
add_library (MinLogLH SHARED ${lib_srcs})
LINK_DIRECTORIES(${MINUIT2}/lib)
target_link_libraries( MinLogLH Minuit2 )

target_include_directories (MinLogLH PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${MINUIT2}/include)