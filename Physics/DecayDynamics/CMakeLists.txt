################################
# Create DecayDynamics library #
################################

set(lib_srcs AbstractDynamicalFunction.cpp RelativisticBreitWigner.cpp
    NonResonant.cpp AmpFlatteRes.cpp Voigtian.cpp Utils/Faddeeva.cc)

set(lib_headers AbstractDynamicalFunction.hpp
    NonResonant.hpp RelativisticBreitWigner.hpp
    AmpFlatteRes.hpp Voigtian.hpp Utils/Faddeeva.hh FormFactor.hpp )

add_library(DecayDynamics
  SHARED ${lib_srcs}
)

target_link_libraries(DecayDynamics
  PUBLIC Core
)

target_include_directories(DecayDynamics
  PUBLIC ${Boost_INCLUDE_DIR}
)

#
# Install
#
install(FILES ${lib_headers}
    DESTINATION include/ComPWA/Physics/DecayDynamics
)
install(TARGETS DecayDynamics
	LIBRARY DESTINATION lib/ComPWA
)