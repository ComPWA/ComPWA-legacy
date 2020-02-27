################################
# Create Dynamics library #
################################

set(lib_srcs RelativisticBreitWigner.cpp Flatte.cpp Voigtian.cpp
             Utils/Faddeeva.cc FormFactor.cpp)

set(lib_headers RelativisticBreitWigner.hpp Flatte.hpp Voigtian.hpp
                Utils/Faddeeva.hh FormFactor.hpp)

add_library(Dynamics ${lib_srcs} ${lib_headers})

target_link_libraries(Dynamics PUBLIC Core FunctionTree)

target_include_directories(Dynamics PUBLIC ${Boost_INCLUDE_DIR})

#
# Install
#
install(FILES ${lib_headers} DESTINATION include/Physics/Dynamics)
install(
  TARGETS Dynamics
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

#
# TESTING
#
#add_subdirectory(test)
