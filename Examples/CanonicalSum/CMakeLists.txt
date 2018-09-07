##########################################################
# Create executable for CanonicalPWA
##########################################################

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
)

target_include_directories( CanonicalPWA
  PUBLIC 
    $<INSTALL_INTERFACE:include>
    ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR}
)

install(TARGETS CanonicalPWA
  RUNTIME DESTINATION bin
)

else ()
  message( WARNING "Required targets not found! Not building\
                    CanonicalPWA executable!")
endif()
