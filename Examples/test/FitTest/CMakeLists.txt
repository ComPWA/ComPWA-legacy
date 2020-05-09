if(TARGET RootData)

  # FitTest
  set(TEST_NAME FitTest)
  add_executable(${TEST_NAME} ${TEST_NAME}.cpp)
  target_include_directories(
    ${TEST_NAME}
    PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
           ${Boost_INCLUDE_DIR}
  )
  target_link_libraries(
    ${TEST_NAME}
    PRIVATE Boost::unit_test_framework
            Core
            Integration
            MinLogLH
            Minuit2IF
            Physics
            RootData
  )
  set_target_properties(
    ${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                            ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/
  )
  add_test(
    NAME ${TEST_NAME}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test
    COMMAND ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/${TEST_NAME}
  )
  configure_file(
    ${PROJECT_SOURCE_DIR}/Physics/particle_list.xml
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/test/particle_list.xml
    COPYONLY
  )

endif()
