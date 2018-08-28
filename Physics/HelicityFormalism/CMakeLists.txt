##############################
# Helicity Formalism library #
##############################

set(lib_srcs 
	../IncoherentIntensity.cpp ../CoherentIntensity.cpp
    ../SequentialPartialAmplitude.cpp AmpWignerD.cpp
    HelicityDecay.cpp HelicityKinematics.cpp
)

set(lib_headers HelicityDecay.hpp AmpWignerD.hpp HelicityKinematics.hpp)

add_library(HelicityFormalism
  SHARED ${lib_srcs}
)

target_link_libraries(HelicityFormalism
  Core
  DecayDynamics
  qft++
)

target_include_directories(HelicityFormalism
	PRIVATE ${QFTPP_INCLUDE_DIR}
)

#
# Install
#
install(FILES ${lib_headers}
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
