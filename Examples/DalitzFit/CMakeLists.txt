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

if( TARGET Minuit2IF AND TARGET MinLogLH
    AND TARGET Tools AND TARGET HelicityFormalism
  )

file( GLOB lib_srcs DalitzFitApp.cpp )
file( GLOB lib_headers ../../Physics/ParticleList.hpp )

add_executable ( DalitzFit
    DalitzFitApp.cpp ${lib_headers} )

target_link_libraries( DalitzFit
    Core 
    Minuit2IF
    MinLogLH
    Tools
    Plotting
    HelicityFormalism
    ${ROOT_LIBRARIES}
    ${Boost_LIBRARIES}
<<<<<<< HEAD
=======
    Threads::Threads #pthread - We need it here, no ideo why!
>>>>>>> master
)

target_include_directories( DalitzFit
  PUBLIC 
    $<INSTALL_INTERFACE:include>
    ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

install(TARGETS DalitzFit
    RUNTIME DESTINATION bin
)

else ()
  message( WARNING "Required targets not found! Not building\
                    DalitzFit executable!")
endif()
