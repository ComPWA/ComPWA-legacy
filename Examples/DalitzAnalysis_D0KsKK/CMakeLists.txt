# Create executable for Dalitz plot analysis of D0->KsKK.

#
# FindRoot
#
find_package( ROOT REQUIRED )

#
# FindBoost
#
find_package(Boost COMPONENTS
               program_options
               REQUIRED QUIET
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
  pthread
)

target_include_directories( DalitzFitD03K
  PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR}
)

install(TARGETS DalitzFitD03K
  RUNTIME DESTINATION bin
)
