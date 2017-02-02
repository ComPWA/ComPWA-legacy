# Create MinLogLH library.
FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)
add_library (MinLogLH SHARED ${lib_srcs} ${lib_headers})
LINK_DIRECTORIES(${MINUIT2}/lib)
target_link_libraries( MinLogLH Minuit2 )

target_include_directories (MinLogLH PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${MINUIT2}/include)