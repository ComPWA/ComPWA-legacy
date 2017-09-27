##########################################################
# Create executable for Dalitz plot analysis of D0->KsKK #
##########################################################

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

add_executable( DalitzFitD03K
  DalitzFitD03K.cpp systematics.hpp
)

target_link_libraries( DalitzFitD03K
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

target_include_directories( DalitzFitD03K
  PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR}
)

install(TARGETS DalitzFitD03K
  RUNTIME DESTINATION bin
)

ELSE ()
  MESSAGE( WARNING "Required targets not found! Not building\
                    DalitzFitD03K executable!")
ENDIF()
