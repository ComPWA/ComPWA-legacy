##########################################################
# Create executable for Dalitz plot analysis of D0->KsKK #
##########################################################

# Make sure all necessary targets are build
if( TARGET Minuit2IF AND
    TARGET MinLogLH AND
    TARGET RootData AND
    TARGET HelicityFormalism
)

add_executable( DalitzFitD03K
  DalitzFitD03K.cpp systematics.hpp
)

target_link_libraries( DalitzFitD03K
  HelicityFormalism
  Minuit2IF
  MinLogLH
  RootData
  Core
  Plotting
  Tools
  ROOT::RIO
  Boost::program_options
)

target_include_directories( DalitzFitD03K
  PUBLIC
    $<INSTALL_INTERFACE:include>
    ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR}
)

install(TARGETS DalitzFitD03K
  RUNTIME DESTINATION bin
)
install(FILES model.xml
  DESTINATION bin
)

else ()
  MESSAGE( WARNING "Required targets not found! Not building\
                    DalitzFitD03K executable!")
endif()
