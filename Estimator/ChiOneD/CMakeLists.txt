# Create ChiOneD library.
SET(lib_srcs ChiOneD.cpp)
SET(lib_headers ChiOneD.hpp ../Estimator.hpp)

add_library(
  ComPWA_ChiOneD
  SHARED ${lib_srcs}  ${lib_headers}
)
target_link_libraries(
  ComPWA_ChiOneD
  ComPWA_Core
)

target_include_directories(
  ComPWA_ChiOneD
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

install (FILES ${lib_headers} DESTINATION include/ComPWA)
