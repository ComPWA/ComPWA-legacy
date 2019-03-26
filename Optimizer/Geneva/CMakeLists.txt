# Geneva minimizer module. Check first if library is available.
if(${GENEVA_FOUND})
  set(lib_srcs GStartIndividual.cpp GenevaIF.cpp GenevaResult.cpp)
  set(lib_headers GStartIndividual.hpp GenevaIF.hpp GenevaResult.hpp)

  add_library(GenevaIF
    SHARED ${lib_srcs} ${lib_headers}
  )
  set_property(TARGET GenevaIF PROPERTY CXX_STANDARD 14)
  
  target_link_libraries(GenevaIF
    PUBLIC Core Boost::serialization
    PRIVATE ${GENEVA_LIBRARIES} Boost::program_options Boost::filesystem
  )

  target_include_directories(GenevaIF
    PUBLIC ${GENEVA_INCLUDE_DIR}
  )

  install(FILES ${lib_headers} DESTINATION include/Optimizer/Geneva)
  install(TARGETS GenevaIF
    LIBRARY DESTINATION lib
  )
endif()
