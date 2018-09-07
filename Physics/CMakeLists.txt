add_subdirectory(DecayDynamics)
add_subdirectory(HelicityFormalism)
add_subdirectory(EvtGen)
add_subdirectory(ExpertSystem)

install(FILES 
	CoherentIntensity.hpp PartialAmplitude.hpp 
    IncoherentIntensity.hpp Amplitude.hpp ParticleList.hpp 
	SequentialPartialAmplitude.hpp
    DESTINATION include/Physics
)
