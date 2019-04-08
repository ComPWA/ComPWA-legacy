################################
# Create Dynamics library #
################################

set(lib_srcs RelativisticBreitWigner.cpp
    NonResonant.cpp FormFactorDecorator.cpp Flatte.cpp Voigtian.cpp Utils/Faddeeva.cc)

set(lib_headers AbstractDynamicalFunction.hpp
    NonResonant.hpp RelativisticBreitWigner.hpp FormFactorDecorator.hpp
    Flatte.hpp Voigtian.hpp Utils/Faddeeva.hh FormFactor.hpp )

add_library(Dynamics
  SHARED ${lib_srcs} ${lib_headers}
)

target_link_libraries(Dynamics
  PUBLIC Core
)

target_include_directories(Dynamics
  PUBLIC ${Boost_INCLUDE_DIR}
)

#
# Install
#
install(FILES ${lib_headers}
  DESTINATION include/Physics/Dynamics
)
install(TARGETS Dynamics
	LIBRARY DESTINATION lib
)

#
# TESTING
#
add_subdirectory(test)
