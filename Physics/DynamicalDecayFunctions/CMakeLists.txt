# Create DynamicalDecayFunctions library.

# The REQUIRED statement is somehow not working
#find_package( ROOT REQUIRED )
#find_package( GSL REQUIRED gsl gslcblas )

FILE(GLOB_RECURSE lib_srcs *.cpp)
FILE(GLOB_RECURSE lib_headers *.hpp)
#FILE(APPEND lib_srcs TwoBodyDecay/RelativisticBreitWigner.cpp)
#FILE(GLOB lib_srcs TwoBodyDecay/TopNodeConstantValue.cpp)
#FILE(GLOB lib_headers TwoBodyDecay/RelativisticBreitWigner.hpp)
#FILE(GLOB lib_headers TwoBodyDecay/TopNodeConstantValue.hpp)

add_library ( ComPWA_DynamicalDecayFunctions SHARED ${lib_srcs} ${lib_headers})
target_link_libraries( ComPWA_DynamicalDecayFunctions ComPWA_DPKinematics ComPWA_Core ${Boost_LIBRARIES} )
target_include_directories( ComPWA_DynamicalDecayFunctions PUBLIC ${Boost_INCLUDE_DIR} )

#install (TARGETS ComPWA_DynamicalDecayFunctions DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)
