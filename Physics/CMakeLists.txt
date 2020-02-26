add_subdirectory(Dynamics)
add_subdirectory(HelicityFormalism)
add_subdirectory(EvtGen)

# Define and install library
set(lib_srcs PhspVolume.cpp)

set(lib_headers PhspVolume.hpp)

add_library(Physics ${lib_srcs} ${lib_headers})

target_link_libraries(Physics PRIVATE Core Integration Data)

install(
  TARGETS Physics
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

# Install headers
install(FILES ${lib_headers} SubSystem.hpp BuilderXML.hpp
              ParticleStateTransitionKinematicsInfo.hpp
        DESTINATION include/Physics)

# Install config files
install(FILES particle_list.xml DESTINATION share)

# Unit tests
add_subdirectory(test)
