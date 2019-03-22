##############################
# Helicity Formalism library #
##############################

set(lib_srcs
	../IncoherentIntensity.cpp ../CoherentIntensity.cpp
    ../SequentialAmplitude.cpp ../SubSystem.cpp
    ../IntensityBuilderXML.cpp ../StrengthIntensityDecorator.cpp
    ../ParticleStateTransitionKinematicsInfo.cpp
    AmpWignerD.cpp HelicityDecay.cpp HelicityKinematics.cpp 
    ../NormalizationIntensityDecorator.cpp
    ../NormalizationAmplitudeDecorator.cpp
)

set(lib_headers AmpWignerD.hpp HelicityDecay.hpp HelicityKinematics.hpp
    ../IncoherentIntensity.hpp ../CoherentIntensity.hpp
    ../SequentialAmplitude.hpp ../SubSystem.hpp
    ../IntensityBuilderXML.hpp ../StrengthIntensityDecorator.hpp
    ../ParticleStateTransitionKinematicsInfo.hpp
    ../NormalizationIntensityDecorator.hpp
    ../NormalizationAmplitudeDecorator.hpp
)

#
# SHARED Library
#
add_library(HelicityFormalism SHARED
 ${lib_srcs} ${lib_headers} 
 )
target_link_libraries(HelicityFormalism
  PUBLIC
    Core
    Data
    Dynamics
    qft++
    Integration
)
target_include_directories(HelicityFormalism
  PRIVATE ${QFTPP_INCLUDE_DIR}
)
install(TARGETS HelicityFormalism
   LIBRARY DESTINATION lib
)

# Install headers
install(FILES ${lib_headers}
    DESTINATION include/Physics/HelicityFormalism
)

#
# TESTING
#
# Testing routines are stored in separate directory
add_subdirectory(test)
