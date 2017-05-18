# Create GenevaIF library.

SET(lib_srcs GStartIndividual.cpp GenevaIF.cpp GenevaResult.cpp)
SET(lib_headers GStartIndividual.hpp GenevaIF.hpp GenevaResult.hpp ../Optimizer.hpp)

add_library( GenevaIF
  SHARED ${lib_srcs} ${lib_headers}
)

target_link_libraries( GenevaIF
  Core
  ${Boost_LIBRARIES}
  ${GENEVA_LIBRARIES}
)

target_include_directories( GenevaIF
  PUBLIC ${GENEVA_INCLUDE_DIR}
)

install (FILES ${lib_headers} DESTINATION include/ComPWA)
install(TARGETS GenevaIF
        LIBRARY DESTINATION lib
)