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
target_compile_definitions(minuitTest PRIVATE "BOOST_TEST_DYN_LINK=1")
add_test(NAME minuitTest
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/minuitTest
)
endif()

if(${GENEVA_FOUND})
add_executable(genevaTest GenevaTestApp.cpp)
target_link_libraries(genevaTest
	Core
    ${GENEVA_LIBRARIES}
)
set_target_properties(genevaTest
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
)
add_test(NAME genevaTest
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/genevaTest
)
endif()

if(${ROOT_MINUIT2_FOUND} AND ${GENEVA_FOUND})
add_executable(optimizerTest OptimizerTestApp.cpp)
target_link_libraries(optimizerTest
	Core
    ROOT::Minuit2
    ${GENEVA_LIBRARIES}
)
set_target_properties(optimizerTest
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
)
add_test(NAME optimizerTest
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/optimizerTest
)
endif()