# Create Minuit2IF library.

# The REQUIRED statement is somehow not working
find_package( ROOT REQUIRED COMPONENTS Minuit2 )
find_package( GSL REQUIRED gsl gslcblas )

FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)
add_library( ComPWA_Minuit2IF SHARED ${lib_srcs} ${lib_headers} ../Optimizer.hpp )
target_link_libraries( ComPWA_Minuit2IF ComPWA_Estimator ComPWA_Core ${Boost_LIBRARIES} ${ROOT_LIBRARIES} ${GSL_LIBRARIES} )
target_include_directories( ComPWA_Minuit2IF PUBLIC ${ROOT_INCLUDE_DIR} ${GSL_INCLUDE_DIR} )

#install (TARGETS ComPWA_Minuit2IF DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)
