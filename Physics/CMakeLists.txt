add_subdirectory( DecayDynamics )
add_subdirectory( HelicityFormalism )
add_subdirectory( qft++ )
#add_subdirectory( EvtGen )

install (FILES CoherentIntensity.hpp PartialAmplitude.hpp 
    IncoherentIntensity.hpp Amplitude.hpp ParticleList.hpp 
	SequentialPartialAmplitude.hpp
    DESTINATION include/ComPWA/Physics
)
