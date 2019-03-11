# Create python module for ComPWA

if(${ROOT_FOUND} AND ${Boost_FOUND})

pybind11_add_module(pycompwa_interface MODULE PyComPWA.cpp)

target_include_directories(pycompwa_interface
  PUBLIC
    $<INSTALL_INTERFACE:include>
    ${ROOT_INCLUDE_DIR} 
    ${PYTHON_INCLUDE_DIRS}
)
target_link_libraries(pycompwa_interface
    PUBLIC
      HelicityFormalism
      Minuit2IF
      MinLogLH
      RootReader
      Tools
      Plotting
      Core
      Boost::serialization
      GSL::gsl
      GSL::gslcblas
      Threads::Threads #pthread 
)

#IF(GENEVA_FOUND)
#  target_link_libraries( PyComPWA
#    GenevaIF
#  )
#  target_include_directories( pycompwa_interface
#    ${GENEVA_INCLUDE_DIR}
#  )
#ENDIF()

#
# Install
#
install(TARGETS pycompwa_interface
    LIBRARY DESTINATION lib
)
endif()
