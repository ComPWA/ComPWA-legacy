# Create executable for Dalitz plot analysis of D0->KsKK.

# The REQUIRED statement is somehow not working
find_package( ROOT REQUIRED )
find_package( Boost REQUIRED )

FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)
add_executable ( DalitzFitD03K ${lib_srcs} ${lib_headers})
target_link_libraries( DalitzFitD03K AmplitudeSum Minuit2 MinLogLH Estimator RootReader DPKinematics Core ${ROOT_LIBRARIES} ${Boost_LIBRARIES} pthread )
target_include_directories( DalitzFitD03K PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

#install (TARGETS DalitzFitD03K DESTINATION bin/ )
#install (FILES ${lib_headers} DESTINATION include/ComPWA/ )
