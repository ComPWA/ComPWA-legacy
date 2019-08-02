add_subdirectory(Dynamics)
add_subdirectory(HelicityFormalism)
add_subdirectory(EvtGen)

install(FILES 
  HelicityFormalism.hpp
  SubSystem.hpp
  ParticleList.hpp
  IntensityBuilderXML.hpp
  ParticleStateTransitionKinematicsInfo.hpp
  DESTINATION include/Physics
)
