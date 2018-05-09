###########################################################
# Create executable for Dalitz plot fit of J/psi -> ypipi #
###########################################################

#
# FindBoost
#
find_package(Boost COMPONENTS
               log
               program_options
               REQUIRED QUIET
)

IF( TARGET Minuit2IF AND TARGET MinLogLH
    AND TARGET Tools AND TARGET HelicityFormalism
  )

FILE( GLOB lib_srcs SimFit.cpp )
FILE( GLOB lib_headers ../../Physics/ParticleList.hpp )

add_executable ( SimFit 
  SimFit.cpp ${lib_headers} )

target_link_libraries( SimFit
    Core 
    Minuit2IF
    MinLogLH
    Tools
    RootReader
    HelicityFormalism
    ${ROOT_LIBRARIES}
    ${Boost_LIBRARIES}
    #pthread # We need it here, no ideo why!
)

target_include_directories( SimFit
    PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

install(TARGETS SimFit
    RUNTIME DESTINATION bin
)

ELSE ()
  MESSAGE( WARNING "Required targets not found! Not building\
                    SimFit executable!")
ENDIF()
