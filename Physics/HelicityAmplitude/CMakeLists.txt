# Create HelicityAmplitude library.

SET(lib_srcs CoherentIntensity.cpp)
SET(lib_headers CoherentIntensity.hpp SequentialTwoBodyDecay.hpp PartialDecay.hpp)

add_library ( ComPWA_HelicityFormalism SHARED ${lib_srcs} ${lib_headers})
target_link_libraries( ComPWA_HelicityFormalism ComPWA_DPKinematics ComPWA_Core ${Boost_LIBRARIES} )
target_include_directories( ComPWA_HelicityFormalism PUBLIC ${Boost_INCLUDE_DIR} )

#install (TARGETS ComPWA_HelicityAmplitude DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)
