add_subdirectory(Dynamics)
add_subdirectory(HelicityFormalism)
add_subdirectory(EvtGen)
# add_subdirectory(ExpertSystem)

install(FILES 
  HelicityFormalism.hpp
	IncoherentIntensity.hpp CoherentIntensity.hpp
  Amplitude.hpp CoefficientAmplitudeDecorator.hpp
	SequentialAmplitude.hpp SubSystem.hpp
	ParticleList.hpp IntensityBuilderXML.hpp
	ParticleStateTransitionKinematicsInfo.hpp
	NormalizationIntensityDecorator.hpp
	NormalizationAmplitudeDecorator.hpp
  StrengthIntensityDecorator.hpp
    DESTINATION include/Physics
)
