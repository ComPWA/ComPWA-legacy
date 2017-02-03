# Create MinLogLH library.

FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)
add_library (MinLogLH SHARED ${lib_srcs} ${lib_headers})

target_link_libraries( MinLogLH Estimator Optimizer )

#target_include_directories (MinLogLH PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${MINUIT2}/include)

install (TARGETS MinLogLH DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)
