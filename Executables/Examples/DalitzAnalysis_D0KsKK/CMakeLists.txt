# Create MinLogLH library.

# The REQUIRED statement is somehow not working
find_package( ROOT REQUIRED )
find_package( GSL REQUIRED gsl gslcblas )

FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)
add_executable ( DalitzFitD03K ${lib_srcs} ${lib_headers})
target_link_libraries( DalitzFitD03K DPKinematics Core ${ROOT_LIBRARIES} ${Boost_LIBRARIES} )
target_include_directories( DalitzFitD03K PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

#install (TARGETS DalitzFitD03K DESTINATION bin/DalitzFitApp )
#install (FILES ${lib_headers} DESTINATION include/ComPWA/ )
