# Create AmplitudeSum library.

# The REQUIRED statement is somehow not working
find_package( ROOT REQUIRED )
find_package( GSL REQUIRED gsl gslcblas )

FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)
add_library ( ComPWA_AmplitudeSum SHARED ${lib_srcs} ${lib_headers})
target_link_libraries( ComPWA_AmplitudeSum ComPWA_Core ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${GSL_LIBRARIES} )
target_include_directories( ComPWA_AmplitudeSum PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} ${GSL_INCLUDE_DIR} )

#install (TARGETS ComPWA_AmplitudeSum DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)
