# Install and export
set(HEADER_LIST HelicityKinematics.hpp WignerD.hpp)
install(FILES ${HEADER_LIST}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${MODULE_NAME}/HelicityFormalism
)

# Add unit tests
add_subdirectory(test)
