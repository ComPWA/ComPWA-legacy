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
    AND TARGET Tools AND TARGET HelicityFormalism
  )

FILE( GLOB lib_srcs Benchmark.cpp )
FILE( GLOB lib_headers ../../Physics/ParticleList.hpp )

add_executable ( Benchmark
    Benchmark.cpp ${lib_headers} )

target_link_libraries( Benchmark
    Core 
    Minuit2IF
    MinLogLH
    Tools
    HelicityFormalism
    ${ROOT_LIBRARIES}
    ${Boost_LIBRARIES}
    #pthread # We need it here, no ideo why!
)

target_include_directories( Benchmark
    PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

install(TARGETS Benchmark
    RUNTIME DESTINATION bin
)

ELSE ()
  MESSAGE( WARNING "Required targets not found! Not building\
                    Benchmark executable!")
ENDIF()
