# Create python module for Dalitz plot analysis

# The REQUIRED statement is somehow not working
find_package( ROOT REQUIRED QUIET )

find_package(Boost COMPONENTS
               program_options
               REQUIRED QUIET
)

#IF(PYTHONLIBS_FOUND)

  #FILE(GLOB lib_srcs *.cpp)
  #FILE(GLOB lib_headers *.hpp)
  add_library ( PythonDalitzFit SHARED 
    PythonFit.cpp PythonFit.hpp
  )
  target_link_libraries( PythonDalitzFit 
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
  IF(GENEVA_FOUND)
      target_link_libraries( PythonDalitzFit 
        GenevaIF
      )
  ENDIF()
  target_include_directories( PythonDalitzFit 
    PUBLIC 
    ${GENEVA_INCLUDE_DIR}
    ${ROOT_INCLUDE_DIR} 
    ${Boost_INCLUDE_DIR}
    ${PYTHON_INCLUDE_DIRS}
  )

  install (TARGETS PythonDalitzFit DESTINATION lib/ )
  install (FILES ${lib_headers} DESTINATION include/ComPWA/ )

  add_subdirectory("${COMPWA_SOURCE_DIR}/Tools/pybind11" "${CMAKE_CURRENT_BINARY_DIR}/pybind11_build")
  pybind11_add_module(Dalitz_ext MODULE Dalitz_ext.cpp)
  target_link_libraries( Dalitz_ext 
    PUBLIC
    Core 
    PythonDalitzFit 
    ${Boost_LIBRARIES}
  )
  target_include_directories( Dalitz_ext 
    PUBLIC 
    ${GENEVA_INCLUDE_DIR}
    ${ROOT_INCLUDE_DIR} 
    ${Boost_INCLUDE_DIR}
    ${PYTHON_INCLUDE_DIRS}
  )

#ENDIF()
