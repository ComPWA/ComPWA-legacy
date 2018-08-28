# Create GenevaIF library.

# Geneva minimizer module. Check first if library is available.
if(${GENEVA_FOUND})
set(lib_srcs GStartIndividual.cpp GenevaIF.cpp GenevaResult.cpp)
set(lib_headers GStartIndividual.hpp GenevaIF.hpp GenevaResult.hpp)

add_library(GenevaIF
  SHARED ${lib_srcs}
)

target_link_libraries(GenevaIF
  PUBLIC Core Boost::serialization
  PRIVATE ${GENEVA_LIBRARIES}
)

target_include_directories(GenevaIF
  PRIVATE ${GENEVA_INCLUDE_DIR}
)

install(FILES ${lib_headers} DESTINATION include/ComPWA)
install(TARGETS GenevaIF
        LIBRARY DESTINATION lib
)
else ()
	message(WARNING "Geneva not found! Not building Geneva interface!")
endif()