if(TARGET Tools)
  add_executable(FormFactorTests FormFactorTests.cpp)

  target_link_libraries(
    FormFactorTests
    Core
    HelicityFormalism
    Tools
    Boost::unit_test_framework
  )

  target_include_directories(
    FormFactorTests PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR}
  )

  # Move testing binaries into a testBin directory
  set_target_properties(
    FormFactorTests
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY
               ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/
  )

  add_test(
    NAME FormFactorTests
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/
    COMMAND
      ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/FormFactorDecoratorTests
  )

else()
  message(WARNING "Requirements not found! Not building tests!")
endif()
