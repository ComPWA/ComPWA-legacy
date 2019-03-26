if(${ROOT_MINUIT2_FOUND})
  add_executable(minuitTest MinuitTestApp.cpp)
  target_link_libraries(minuitTest
	Core
    Minuit2IF
    Boost::unit_test_framework
  )
  set_target_properties(minuitTest
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
  )
  add_test(NAME minuitTest
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/minuitTest
  )
endif()

if(${GENEVA_FOUND})
  add_executable(genevaTest GenevaTestApp.cpp)
  target_link_libraries(genevaTest
    Core
	GenevaIF
	Boost::unit_test_framework
  )
  target_include_directories(genevaTest
    PUBLIC ${GENEVA_INCLUDE_DIR}
  )
  set_target_properties(genevaTest
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
  )
  add_test(NAME genevaTest
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/genevaTest
  )
endif()