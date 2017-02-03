# Create Optimizer library.
FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)
#add_library (Physics SHARED ${lib_srcs} ${lib_headers})

add_subdirectory( AmplitudeSum )
add_subdirectory( DPKinematics )

#target_link_libraries( Physics Core ${Boost_LIBRARIES} )

#target_include_directories( Physics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

#install (TARGETS Physics DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)
