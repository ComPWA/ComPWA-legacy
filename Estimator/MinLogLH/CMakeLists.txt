# Create MinLogLH library.

SET(lib_srcs MinLogLH.cpp)
SET(lib_headers MinLogLH.hpp ../Estimator.hpp)

add_library(
  ComPWA_MinLogLH
  SHARED ${lib_srcs} ${lib_headers}
)

target_link_libraries(
  ComPWA_MinLogLH
  ComPWA_Core
)

install (FILES ${lib_headers} DESTINATION include/ComPWA)
