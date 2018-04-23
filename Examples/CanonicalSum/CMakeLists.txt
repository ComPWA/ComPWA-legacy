##########################################################
# Create executable for CanonicalPWA
##########################################################

#
# FindBoost
#
find_package(Boost COMPONENTS
               program_options
               REQUIRED QUIET
)

# Make sure all necessary targets are build
if( TARGET Minuit2IF AND TARGET MinLogLH AND TARGET RootReader
    AND TARGET Tools AND TARGET HelicityFormalism
)

add_executable( CanonicalPWA
  CanonicalPWA.cpp
)

target_link_libraries( CanonicalPWA
  HelicityFormalism
  Minuit2IF
  MinLogLH
  RootReader
  Tools
  Core
  ${ROOT_LIBRARIES}
  ${Boost_LIBRARIES}
  pthread # We need it here, no ideo why!
)

target_include_directories( CanonicalPWA
  PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR}
)

install(TARGETS CanonicalPWA
  RUNTIME DESTINATION bin
)

ELSE ()
  MESSAGE( WARNING "Required targets not found! Not building\
                    CanonicalPWA executable!")
ENDIF()
