add_executable(Optimizer_MinuitTest MinuitTestApp.cpp)
target_link_libraries(Optimizer_MinuitTest Core Minuit2IF
                      Boost::unit_test_framework)
set_target_properties(
  Optimizer_MinuitTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                  ${PROJECT_BINARY_DIR}/bin/test/)
add_test(
  NAME Optimizer_MinuitTest
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
  COMMAND ${PROJECT_BINARY_DIR}/bin/test/Optimizer_MinuitTest)

add_executable(Optimizer_MinuitFitResultSerializationTest
               FitResultSerializationTest.cpp)
target_link_libraries(Optimizer_MinuitFitResultSerializationTest Core Minuit2IF
                      Boost::unit_test_framework)
set_target_properties(
  Optimizer_MinuitFitResultSerializationTest
  PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/)
add_test(
  NAME Optimizer_MinuitFitResultSerializationTest
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
  COMMAND
    ${PROJECT_BINARY_DIR}/bin/test/Optimizer_MinuitFitResultSerializationTest)
