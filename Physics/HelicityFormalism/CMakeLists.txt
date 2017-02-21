# Create HelicityAmplitude library.

SET( lib_srcs IncoherentIntensity.cpp CoherentIntensity.cpp AbstractDynamicalFunction.cpp RelativisticBreitWigner.cpp AmpWignerD.cpp PartialDecay.cpp )
SET( lib_headers IncoherentIntensity.hpp CoherentIntensity.hpp SequentialTwoBodyDecay.hpp PartialDecay.hpp AbstractDynamicalFunction.hpp NonResonant.hpp PartialDecay.hpp RelativisticBreitWigner.hpp SequentialTwoBodyDecay.hpp AmpWignerD.hpp )

add_library ( ComPWA_HelicityFormalism SHARED ${lib_srcs} ${lib_headers})
target_link_libraries( ComPWA_HelicityFormalism ComPWA_Core ${Boost_LIBRARIES} )
target_include_directories( ComPWA_HelicityFormalism PUBLIC ${Boost_INCLUDE_DIR} )

#install (TARGETS ComPWA_HelicityAmplitude DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)
