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

IF( TARGET Minuit2IF AND TARGET MinLogLH
    AND TARGET Tools AND TARGET EvtGenIF
  )

FILE( GLOB lib_srcs EvtGen.cpp )
FILE( GLOB lib_headers ../../Physics/ParticleList.hpp )

add_executable ( EvtGen
    EvtGen.cpp ${lib_headers} )

target_link_libraries( EvtGen
    Core 
    Minuit2IF
    MinLogLH
    Tools
    EvtGenIF
    ${ROOT_LIBRARIES}
    ${Boost_LIBRARIES}
    pthread # We need it here, no ideo why!
)

target_include_directories( EvtGen
    PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

install(TARGETS EvtGen
    RUNTIME DESTINATION bin
)

ELSE ()
  MESSAGE( WARNING "Required targets not found! Not building\
                    EvtGen executable!")
ENDIF()
