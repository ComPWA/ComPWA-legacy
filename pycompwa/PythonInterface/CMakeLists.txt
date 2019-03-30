# Create python module for ComPWA

if(${ROOT_FOUND} AND ${Boost_FOUND})

pybind11_add_module(ui MODULE PyComPWA.cpp)

target_include_directories(ui
  PUBLIC
    $<INSTALL_INTERFACE:include>
    ${ROOT_INCLUDE_DIR} 
    ${PYTHON_INCLUDE_DIRS}
)

target_link_libraries(ui
  PRIVATE Core Data RootDataIO MinLogLH Minuit2IF HelicityFormalism Tools
  Plotting Boost::serialization Threads::Threads
)

add_custom_target(copy ALL DEPENDS ui)

add_custom_command(TARGET copy
  COMMAND ${CMAKE_COMMAND} -E copy
  $<TARGET_FILE:ui> ${CMAKE_CURRENT_BINARY_DIR}/../pycompwa/
)

#
# Install
#
install(TARGETS ui
    LIBRARY DESTINATION lib
)
endif()
