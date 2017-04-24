# Create executable for Dalitz plot analysis of D0->KsKK.

# The REQUIRED statement is somehow not working
find_package( ROOT REQUIRED )
find_package( Boost REQUIRED )

SET(lib_srcs DalitzFitD03K.cpp)
SET(lib_headers Tools.hpp)

add_executable (
  DalitzFitD03K ${lib_srcs} ${lib_headers}
)

target_link_libraries(
  DalitzFitD03K
  ComPWA_HelicityFormalism
  ComPWA_Minuit2IF
  ComPWA_MinLogLH
  ComPWA_RootReader
  ComPWA_Core
  ${ROOT_LIBRARIES}
  ${Boost_LIBRARIES}
  pthread
)

target_include_directories(
  DalitzFitD03K
  PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR}
)

#install (TARGETS DalitzFitD03K DESTINATION bin/ )
#install (FILES ${lib_headers} DESTINATION include/ComPWA/ )
