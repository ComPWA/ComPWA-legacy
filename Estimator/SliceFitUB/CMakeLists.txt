# Create SliceFitUB library.
find_package( ROOT REQUIRED )

SET(lib_srcs SliceFitUB.cpp)
SET(lib_headers SliceFitUB.hpp ../Estimator.hpp)

add_library(
  ComPWA_SliceFitUB
  SHARED ${lib_srcs} ${lib_headers}
)

target_link_libraries(
  ComPWA_SliceFitUB
  ComPWA_Core
  ${ROOT_LIBRARIES}
)

target_include_directories(
  ComPWA_SliceFitUB
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${MINUIT2}/include ${ROOT_INCLUDE_DIR}
)

install (FILES ${lib_headers} DESTINATION include/ComPWA)
