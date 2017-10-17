add_subdirectory( DecayDynamics )
add_subdirectory( HelicityFormalism )
add_subdirectory( qft++ )
#add_subdirectory( EvtGen )

install (FILES Resonance.hpp Amplitude.hpp ParticleList.hpp
    DESTINATION include/ComPWA/Physics
)
