# Create HelicityAmplitude library.

FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)
add_library ( ComPWA_HelicityAmplitude SHARED ${lib_srcs} ${lib_headers})
target_link_libraries( ComPWA_HelicityAmplitude ComPWA_DPKinematics ComPWA_Core ${Boost_LIBRARIES} )
target_include_directories( ComPWA_HelicityAmplitude PUBLIC ${Boost_INCLUDE_DIR} )

#install (TARGETS ComPWA_HelicityAmplitude DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)
