# Create python module for Dalitz plot analysis

# The REQUIRED statement is somehow not working
find_package( ROOT REQUIRED )
find_package( Boost REQUIRED )

#FILE(GLOB lib_srcs *.cpp)
#FILE(GLOB lib_headers *.hpp)
add_library ( PythonDalitzFit SHARED 
  PythonFit.cpp PlotData.cpp PythonFit.hpp PlotData.hpp Tools.hpp systematics.hpp
)
target_link_libraries( PythonDalitzFit 
  HelicityFormalism
  Minuit2IF
  MinLogLH
  RootReader
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
)

install (TARGETS PythonDalitzFit DESTINATION lib/ )
install (FILES ${lib_headers} DESTINATION include/ComPWA/ )

IF(PYTHONLIBS_FOUND)
    # Define the wrapper library that wraps our library
    add_library( Dalitz_ext SHARED Dalitz_ext.cpp )
    target_link_libraries( Dalitz_ext Core PythonDalitzFit ${Boost_LIBRARIES} )
    # don't prepend wrapper library name with lib
    set_target_properties( Dalitz_ext PROPERTIES PREFIX "" )
ENDIF()