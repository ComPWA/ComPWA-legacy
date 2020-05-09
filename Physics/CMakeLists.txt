add_subdirectory(Dynamics)
add_subdirectory(HelicityFormalism)
add_subdirectory(EvtGen)

# Define and install library
set(SOURCE_LIST PhspVolume.cpp)

set(HEADER_LIST PhspVolume.hpp)

add_library(Physics ${SOURCE_LIST} ${HEADER_LIST})

target_link_libraries(
  Physics
  PRIVATE Core
          Integration
          Data
)

install(
  TARGETS Physics
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# Install headers
install(
  FILES ${HEADER_LIST}
        SubSystem.hpp
        BuilderXML.hpp
        ParticleStateTransitionKinematicsInfo.hpp
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Physics
)

# Install config files
install(FILES particle_list.xml DESTINATION share)

# Unit tests
add_subdirectory(test)
