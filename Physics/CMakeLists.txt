# Create Physics libraries.
#FILE(GLOB lib_srcs *.cpp)
#FILE(GLOB lib_headers *.hpp)
#add_library (ComPWA_Physics SHARED ${lib_srcs} ${lib_headers})

add_subdirectory( AmplitudeSum )
#add_subdirectory( BreitWigner )
#add_subdirectory( DecayTree )
add_subdirectory( HelicityFormalism )

#target_link_libraries( ComPWA_Physics ComPWA_Core ${Boost_LIBRARIES} )

#target_include_directories( ComPWA_Physics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

#install (TARGETS Physics DESTINATION lib/ComPWA)
#install (FILES ${lib_headers} DESTINATION include/ComPWA)
