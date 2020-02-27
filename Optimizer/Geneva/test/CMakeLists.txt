add_executable(Optimizer_GenevaTest GenevaTestApp.cpp)
target_link_libraries(Optimizer_GenevaTest Core GenevaIF
                      Boost::unit_test_framework)
target_include_directories(Optimizer_GenevaTest PUBLIC ${GENEVA_INCLUDE_DIR})
set_target_properties(
  Optimizer_GenevaTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                  ${PROJECT_BINARY_DIR}/bin/test/)
add_test(
  NAME Optimizer_GenevaTest
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
  COMMAND ${PROJECT_BINARY_DIR}/bin/test/Optimizer_GenevaTest)
