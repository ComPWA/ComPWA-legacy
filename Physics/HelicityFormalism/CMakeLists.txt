##############################
# Helicity Formalism library #
##############################

#
# Find Boost
#
find_package( Boost COMPONENTS
               unit_test_framework
               system
               REQUIRED QUIET
)

SET( lib_srcs ../IncoherentIntensity.cpp ../CoherentIntensity.cpp
    ../SequentialPartialAmplitude.cpp AmpWignerD.cpp
    HelicityDecay.cpp HelicityFromCanonicalSum.cpp HelicityKinematics.cpp )

SET( lib_headers HelicityDecay.hpp AmpWignerD.hpp HelicityKinematics.hpp)

SET( lib_physics_headers ../IncoherentIntensity.hpp ../CoherentIntensity.hpp
    ../SequentialPartialAmplitude.hpp ../Amplitude.hpp ../PartialAmplitude.hpp
    ../ParticleList.hpp  HelicityDecay.hpp 
    AmpWignerD.hpp HelicityFromCanonicalSum.hpp HelicityKinematics.hpp )

add_library( HelicityFormalism
  SHARED ${lib_srcs} ${lib_headers} ${lib_physics_headers}
)

target_link_libraries( HelicityFormalism
  Core
  qft++
  DecayDynamics
  ${Boost_LIBRARIES}
)

target_include_directories( HelicityFormalism
  PUBLIC ${Boost_INCLUDE_DIR}
)

#
# Install
#
install (FILES ${lib_headers}
    DESTINATION include/ComPWA/Physics/HelicityFormalism
)

install(TARGETS HelicityFormalism
    LIBRARY DESTINATION lib/ComPWA
)

#
# TESTING
#
# Testing routines are stored in separate directory
add_subdirectory(test)
