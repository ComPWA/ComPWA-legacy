# Create executable for Dalitz plot analysis of D0->KsKK.

# The REQUIRED statement is somehow not working
find_package( ROOT REQUIRED )
find_package( Boost REQUIRED )

add_library( DalitzFitD03K
  PlotData.cpp PlotData.hpp
)

target_link_libraries( ComPWA_DalitzFitD03K
  ComPWA_HelicityFormalism
  ComPWA_Minuit2IF
  ComPWA_MinLogLH
  ComPWA_RootReader
  ComPWA_Core
  ${ROOT_LIBRARIES}
  ${Boost_LIBRARIES}
  pthread
)

add_executable( DalitzFitD03KApp
  DalitzFitD03K.cpp Tools.hpp
)

target_link_libraries( DalitzFitD03KApp
  ComPWA_DalitzFitD03K
  ComPWA_HelicityFormalism
  ComPWA_Minuit2IF
  ComPWA_MinLogLH
  ComPWA_RootReader
  ComPWA_Core
  ${ROOT_LIBRARIES}
  ${Boost_LIBRARIES}
  pthread
)

target_include_directories( DalitzFitD03K
  PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR}
)
