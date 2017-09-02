# Create executable for Dalitz plot analysis of D0->KsKK.

# The REQUIRED statement is somehow not working
find_package( ROOT REQUIRED )
find_package( Boost REQUIRED )

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
