# Create DecayTree library.

# The REQUIRED statement is somehow not working
#find_package( ROOT REQUIRED )
#find_package( GSL REQUIRED gsl gslcblas )

FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)
add_library ( ComPWA_DecayTree SHARED ${lib_srcs} ${lib_headers})
target_link_libraries( ComPWA_DecayTree ComPWA_DPKinematics ComPWA_Core ${Boost_LIBRARIES} )
target_include_directories( ComPWA_DecayTree PUBLIC ${Boost_INCLUDE_DIR} )

#install (TARGETS ComPWA_DecayTree DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)
