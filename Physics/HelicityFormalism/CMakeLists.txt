##############################
# Helicity Formalism library #
##############################

# Define and install library
set(SOURCE_LIST WignerD.cpp HelicityKinematics.cpp ../SubSystem.cpp
                ../BuilderXML.cpp ../ParticleStateTransitionKinematicsInfo.cpp
)

set(HEADER_LIST WignerD.hpp HelicityKinematics.hpp ../SubSystem.hpp
                ../BuilderXML.hpp ../ParticleStateTransitionKinematicsInfo.hpp
)

add_library(HelicityFormalism ${SOURCE_LIST} ${HEADER_LIST})

target_link_libraries(
  HelicityFormalism
  PUBLIC Core FunctionTree Data Dynamics
  PRIVATE qft++ Integration
)

install(
  TARGETS HelicityFormalism
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

# Install headers
install(FILES ${HEADER_LIST} DESTINATION include/Physics/HelicityFormalism)

# Unit tests
add_subdirectory(test)
