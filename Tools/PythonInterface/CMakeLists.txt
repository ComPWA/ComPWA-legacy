# Create python module for ComPWA

if(${ROOT_FOUND} AND ${Boost_FOUND})

pybind11_add_module(pycompwa MODULE PyComPWA.cpp)

target_include_directories(pycompwa
  PUBLIC
    $<INSTALL_INTERFACE:include>
    ${ROOT_INCLUDE_DIR} 
    ${PYTHON_INCLUDE_DIRS}
)

target_link_libraries(pycompwa
  PRIVATE Core Data RootDataIO MinLogLH Minuit2IF HelicityFormalism Tools
  Plotting Boost::serialization Threads::Threads
)

#
# Install
#
install(TARGETS pycompwa
    LIBRARY DESTINATION lib
)
endif()
