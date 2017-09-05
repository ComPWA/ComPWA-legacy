# Create SliceFit library.
find_package( ROOT REQUIRED QUIET )

SET(lib_srcs SliceFit.cpp)
SET(lib_headers SliceFit.hpp ../Estimator.hpp)

add_library(
  ComPWA_SliceFit
  SHARED ${lib_srcs} ${lib_headers}
)

target_link_libraries(
  ComPWA_SliceFit
  ComPWA_Core
  ${ROOT_LIBRARIES}
)

target_include_directories(
  ComPWA_SliceFit
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${MINUIT2}/include ${ROOT_INCLUDE_DIR}
)

install (FILES ${lib_headers} DESTINATION include/ComPWA)
