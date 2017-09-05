###########################################################
# Create executable for Dalitz plot fit of J/psi -> ypipi #
###########################################################

#
# FindRoot
#
find_package( ROOT REQUIRED QUIET )

#
# FindBoost
#
find_package(Boost COMPONENTS
               program_options
               REQUIRED QUIET
)

FILE( GLOB lib_srcs DalitzFitApp.cpp )
FILE( GLOB lib_headers ../../Physics/ParticleList.hpp )

add_executable ( DalitzFit
    DalitzFitApp.cpp ${lib_headers} )

target_link_libraries( DalitzFit
    Core 
    Minuit2IF
    MinLogLH 
    RootReader
    Tools
    HelicityFormalism
    ${ROOT_LIBRARIES}
    ${Boost_LIBRARIES} 
)

target_include_directories( DalitzFit
    PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )
