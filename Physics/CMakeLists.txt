set(MODULE_NAME Physics)

# Define library
set(SOURCE_LIST
    BuilderXML.cpp
    HelicityFormalism/HelicityKinematics.cpp
    HelicityFormalism/WignerD.cpp
    ParticleStateTransitionKinematicsInfo.cpp
    PhspVolume.cpp
    SubSystem.cpp
)
add_library(${MODULE_NAME} ${SOURCE_LIST})
target_include_directories(
  ${MODULE_NAME}
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
         ${QFTPP_INCLUDE_DIR}
)
target_link_libraries(
  ${MODULE_NAME}
  PRIVATE Core
          Data
          Dynamics
          Integration
)

# Install and export
set(HEADER_LIST
    BuilderXML.hpp
    ParticleStateTransitionKinematicsInfo.hpp
    PhspVolume.hpp
    SubSystem.hpp
)
install(FILES ${HEADER_LIST}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${MODULE_NAME}
)
install(
  TARGETS ${MODULE_NAME}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(FILES particle_list.xml DESTINATION ${CMAKE_INSTALL_DATADIR})

# Add unit tests
# add_subdirectory(test)

# Add subdirectories
add_subdirectory(Dynamics)
add_subdirectory(EvtGen)
add_subdirectory(HelicityFormalism)
