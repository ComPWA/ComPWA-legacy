if(GenevaIF)

  # GenevaTestApp
  set(TEST_NAME GenevaTestApp)
  add_executable(${TEST_NAME} ${TEST_NAME}.cpp)
  target_include_directories(
    ${TEST_NAME}
    PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
           ${GENEVA_INCLUDE_DIR}
  )
  target_link_libraries(
    ${TEST_NAME}
    PRIVATE Core
            GenevaIF
            Boost::unit_test_framework
  )
  set_target_properties(
    ${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                            ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/
  )
  add_test(
    NAME ${TEST_NAME}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/
    COMMAND ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/${TEST_NAME}
  )

endif()
