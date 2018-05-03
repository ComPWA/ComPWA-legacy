# Create python module for ComPWA

find_package( ROOT REQUIRED QUIET )

find_package( Boost COMPONENTS
               program_options
               REQUIRED QUIET
)

# Make sure pybind uses c++11
set(PYBIND11_CPP_STANDARD -std=c++11)
add_subdirectory( pybind11 )

#add_subdirectory( "${COMPWA_SOURCE_DIR}/Tools/pybind11"
#    "${CMAKE_CURRENT_BINARY_DIR}/pybind11_build"
#)

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
