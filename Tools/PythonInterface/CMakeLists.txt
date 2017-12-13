# Create python module for Dalitz plot analysis

# The REQUIRED statement is somehow not working
find_package( ROOT REQUIRED QUIET )

find_package(Boost COMPONENTS
               program_options
               REQUIRED QUIET
)

#IF(PYTHONLIBS_FOUND)

  add_subdirectory("${COMPWA_SOURCE_DIR}/Tools/pybind11" "${CMAKE_CURRENT_BINARY_DIR}/pybind11_build")
  pybind11_add_module(PyComPWA MODULE PyComPWA.cpp)
  target_link_libraries( PyComPWA 
    PUBLIC
    HelicityFormalism
    Minuit2IF
    MinLogLH
    RootReader
    Tools
    Core
    ${ROOT_LIBRARIES} 
    ${Boost_LIBRARIES} 
    pthread 
  )
  target_include_directories( PyComPWA 
    PUBLIC 
    ${ROOT_INCLUDE_DIR} 
    ${Boost_INCLUDE_DIR}
    ${PYTHON_INCLUDE_DIRS}
  )
  #IF(GENEVA_FOUND)
  #  target_link_libraries( PyComPWA 
  #    GenevaIF
  #  )
  #  target_include_directories( PyComPWA  
  #    ${GENEVA_INCLUDE_DIR}
  #  )
  #ENDIF()

#ENDIF()
