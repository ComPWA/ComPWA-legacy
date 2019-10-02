add_subdirectory(Dynamics)
add_subdirectory(HelicityFormalism)
add_subdirectory(EvtGen)

install(FILES 
  HelicityFormalism.hpp
  SubSystem.hpp
  BuilderXML.hpp
  ParticleStateTransitionKinematicsInfo.hpp
  DESTINATION include/Physics
)

install(FILES 
  particle_list.xml
  DESTINATION share
)