##############################
# Helicity Formalism library #
##############################

set(lib_srcs
  WignerD.cpp
  HelicityKinematics.cpp
  ../SubSystem.cpp
  ../IntensityBuilderXML.cpp
  ../ParticleStateTransitionKinematicsInfo.cpp
)

set(lib_headers
  WignerD.hpp
  HelicityKinematics.hpp
  ../SubSystem.hpp
  ../IntensityBuilderXML.hpp
  ../ParticleStateTransitionKinematicsInfo.hpp
)

#
# SHARED Library
#
add_library(HelicityFormalism SHARED
  ${lib_srcs} ${lib_headers} 
)

target_link_libraries(HelicityFormalism
  PUBLIC Core FunctionTree Data Dynamics
  PRIVATE qft++ Integration
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
