# Create python module for ComPWA

if(${ROOT_FOUND} AND ${Boost_FOUND})

pybind11_add_module(pycompwa MODULE PyComPWA.cpp)

target_include_directories(pycompwa
    PUBLIC
    	${ROOT_INCLUDE_DIR} 
    	${Boost_INCLUDE_DIR}
    	${PYTHON_INCLUDE_DIRS}
)
target_link_libraries(pycompwa
    PUBLIC
      HelicityFormalism
      Minuit2IF
      MinLogLH
      RootReader
      Tools
      Plotting
      Core
      Boost::filesystem
      Boost::serialization
      pthread 
)

#IF(GENEVA_FOUND)
#  target_link_libraries( PyComPWA
#    GenevaIF
#  )
#  target_include_directories( pycompwa
#    ${GENEVA_INCLUDE_DIR}
#  )
#ENDIF()

#
# Install
#
install(TARGETS pycompwa
    LIBRARY DESTINATION lib/ComPWA
)
endif()
