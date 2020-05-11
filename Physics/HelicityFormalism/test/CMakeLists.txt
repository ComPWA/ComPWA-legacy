if(TARGET Tools)
  # -------------------- Model Function Tree Tests -------------------- #
  add_executable(ModelFunctionTreeTests ModelFunctionTreeTests.cpp)

  target_link_libraries(
    ModelFunctionTreeTests
    Boost::unit_test_framework
    Core
    HelicityFormalism
    RootData
  )

  target_include_directories(
    ModelFunctionTreeTests PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR}
  )

  # Move testing binaries into a testBin directory
  set_target_properties(
    ModelFunctionTreeTests PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                      ${PROJECT_BINARY_DIR}/bin/test/
  )

  add_test(
    NAME ModelFunctionTreeTests
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/ModelFunctionTreeTests
  )

  # -------------------- Helicity Angle Test -------------------- #
  add_executable(HelicityAngle HelicityAngle.cpp)

  target_link_libraries(
    HelicityAngle
    Boost::unit_test_framework
    Core
    HelicityFormalism
    RootData
    qft++
  )

  target_include_directories(
    HelicityAngle
    PUBLIC ${ROOT_INCLUDE_DIR}
           ${Boost_INCLUDE_DIR}
           ${QFTPP_INCLUDE_DIR}
  )

  # Move testing binaries into a testBin directory
  set_target_properties(
    HelicityAngle PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                             ${PROJECT_BINARY_DIR}/bin/test/
  )

  add_test(
    NAME HelicityAngle
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/HelicityAngle
  )

  # -------------------- Helicity Angles Correctness Test -------------------- #
  add_executable(
    HelicityAnglesCorrectnessTest HelicityAnglesCorrectnessTest.cpp
  )

  target_link_libraries(
    HelicityAnglesCorrectnessTest
    Boost::unit_test_framework
    Core
    HelicityFormalism
    RootData
    qft++
  )

  target_include_directories(
    HelicityAnglesCorrectnessTest
    PUBLIC ${ROOT_INCLUDE_DIR}
           ${Boost_INCLUDE_DIR}
           ${QFTPP_INCLUDE_DIR}
  )

  # Move testing binaries into a testBin directory
  set_target_properties(
    HelicityAnglesCorrectnessTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                             ${PROJECT_BINARY_DIR}/bin/test/
  )

  add_test(
    NAME HelicityAnglesCorrectnessTest
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/HelicityAnglesCorrectnessTest
  )

  # -------------------- Helicity Kinematics Tests -------------------- #
  add_executable(HelicityKinematicsTests HelicityKinematicsTests.cpp)

  target_link_libraries(
    HelicityKinematicsTests
    Boost::unit_test_framework
    Core
    HelicityFormalism
  )

  target_include_directories(
    HelicityKinematicsTests PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR}
  )

  # Move testing binaries into a testBin directory
  set_target_properties(
    HelicityKinematicsTests PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                       ${PROJECT_BINARY_DIR}/bin/test/
  )

  add_test(
    NAME HelicityKinematicsTests
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/HelicityKinematicsTests
  )
else()
  message(WARNING "Requirements not found! Not building tests!")
endif()
