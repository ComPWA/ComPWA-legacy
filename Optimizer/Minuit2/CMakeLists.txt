# Create MinLogLH library.

# The REQUIRED statement is somehow not working
find_package( ROOT REQUIRED COMPONENTS Minuit2 )
find_package( GSL REQUIRED gsl gslcblas )

FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)
add_library (Minuit2 SHARED ${lib_srcs} ${lib_headers})
target_link_libraries( Minuit2 Estimator Optimizer ${ROOT_LIBRARIES} ${GSL_LIBRARIES} )
target_include_directories( Minuit2 PUBLIC ${ROOT_INCLUDE_DIR} ${GSL_INCLUDE_DIR} )

install (TARGETS Minuit2 DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)
