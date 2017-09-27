###########################################################
# Create executable for Dalitz plot fit of J/psi -> ypipi #
###########################################################

#
# FindBoost
#
find_package(Boost COMPONENTS
               program_options
               REQUIRED QUIET
)

IF( TARGET Minuit2IF AND TARGET MinLogLH AND TARGET RootReader
    AND TARGET Tools AND TARGET HelicityFormalism
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
    pthread # We need it here, no ideo why!
)

target_include_directories( DalitzFit
    PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

ELSE ()
  MESSAGE( WARNING "Required targets not found! Not building\
                    DalitzFit executable!")
ENDIF()
