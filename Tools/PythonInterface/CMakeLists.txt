# Create python module for ComPWA

find_package( ROOT REQUIRED QUIET )

find_package( Boost REQUIRED QUIET )

pybind11_add_module( pycompwa MODULE PyComPWA.cpp )

target_link_libraries( pycompwa
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

target_include_directories( pycompwa
    PUBLIC 
    ${ROOT_INCLUDE_DIR} 
    ${Boost_INCLUDE_DIR}
    ${PYTHON_INCLUDE_DIRS}
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
install( TARGETS pycompwa
    LIBRARY DESTINATION lib/ComPWA
)
