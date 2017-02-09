# Create MinLogLH library.

FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)
add_library( ComPWA_MinLogLH SHARED ${lib_srcs} ${lib_headers} )

target_link_libraries( ComPWA_MinLogLH ComPWA_Estimator ComPWA_Core )

#target_include_directories (ComPWA_MinLogLH PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${MINUIT2}/include)

#install (TARGETS ComPWA_MinLogLH DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)
