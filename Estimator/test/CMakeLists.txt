add_executable(Estimator_MinLogLHEstimatorTest MinLogLHEstimatorTest.cpp)

target_link_libraries(
  Estimator_MinLogLHEstimatorTest
  PUBLIC Boost::unit_test_framework
         Core
         Integration
         MinLogLH
         Minuit2IF
         Tools
)

set_target_properties(
  Estimator_MinLogLHEstimatorTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                             ${PROJECT_BINARY_DIR}/bin/test/
)

add_test(
  NAME Estimator_MinLogLHEstimatorTest
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
  COMMAND ${PROJECT_BINARY_DIR}/bin/test/Estimator_MinLogLHEstimatorTest
)
