##############################
# Helicity Formalism library #
##############################

# Define and install library
set(lib_srcs WignerD.cpp HelicityKinematics.cpp ../SubSystem.cpp
             ../BuilderXML.cpp ../ParticleStateTransitionKinematicsInfo.cpp)

set(lib_headers WignerD.hpp HelicityKinematics.hpp ../SubSystem.hpp
                ../BuilderXML.hpp ../ParticleStateTransitionKinematicsInfo.hpp)

add_library(HelicityFormalism ${lib_srcs} ${lib_headers})

target_link_libraries(
  HelicityFormalism
  PUBLIC Core FunctionTree Data Dynamics
  PRIVATE qft++ Integration)

install(
  TARGETS HelicityFormalism
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

# Install headers
install(FILES ${lib_headers} DESTINATION include/Physics/HelicityFormalism)

# Unit tests
add_subdirectory(test)
