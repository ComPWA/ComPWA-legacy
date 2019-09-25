set(testName "CoreTest_ParticleTest")

add_executable(${testName} ParticleTest.cpp)

target_link_libraries(${testName}
  Core
  Threads::Threads
  Boost::unit_test_framework
)

target_include_directories(${testName}
  PRIVATE Boost::unit_test_framework)

# Move testing binaries into a testBin directory
set_target_properties(${testName}
  PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
)
  
add_test(NAME ${testName}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/${testName}
)