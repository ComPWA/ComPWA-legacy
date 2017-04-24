# Create Minuit2IF library.

# The REQUIRED statement is somehow not working
find_package( ROOT REQUIRED COMPONENTS Minuit2 )
find_package( GSL REQUIRED gsl gslcblas )

SET(lib_srcs MinuitFcn.cpp MinuitIF.cpp MinuitResult.cpp)
SET(lib_headers MinuitFcn.hpp MinuitIF.hpp MinuitResult.hpp ../Optimizer.hpp)

add_library( ComPWA_Minuit2IF
  SHARED ${lib_srcs} ${lib_headers}
)

target_link_libraries( ComPWA_Minuit2IF
  ComPWA_Core
  ${Boost_LIBRARIES}
  ${ROOT_LIBRARIES}
  ${GSL_LIBRARIES}
)

target_include_directories( ComPWA_Minuit2IF
  PUBLIC ${ROOT_INCLUDE_DIR} ${GSL_INCLUDE_DIR}
)

install (FILES ${lib_headers} DESTINATION include/ComPWA)
