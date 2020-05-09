set(MODULE_NAME Estimator)

# Install and export
set(HEADER_LIST Estimator.hpp)
install(FILES ${HEADER_LIST}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${MODULE_NAME}
)

# Add subdirectories
add_subdirectory(ChiOneD)
add_subdirectory(MinLogLH)

# Add unit tests
add_subdirectory(test)
