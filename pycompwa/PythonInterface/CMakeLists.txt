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
  PRIVATE Core FunctionTree Data RootDataIO MinLogLH Minuit2IF HelicityFormalism Tools
  Plotting Boost::serialization Threads::Threads ROOT::Hist
)

add_custom_command(TARGET ui
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy
  $<TARGET_FILE:ui> ${CMAKE_CURRENT_BINARY_DIR}/../pycompwa/
)

install(TARGETS ui
  LIBRARY DESTINATION share/ComPWA/pycompwa/pycompwa
)
endif()
