if(TARGET Physics)
  # ------- PhspVolume Tests MC integration ------- #
  add_executable(PhspVolumeTest PhspVolumeTest.cpp)
  target_link_libraries(PhspVolumeTest Boost::unit_test_framework Physics)
  target_include_directories(PhspVolumeTest PUBLIC ${Boost_INCLUDE_DIR})
  set_target_properties(
    PhspVolumeTest PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                              ${PROJECT_BINARY_DIR}/bin/test/)
  add_test(
    NAME PhspVolumeTest
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/PhspVolumeTest)

else()
  message(WARNING "Requirements not found! Not building tests!")
endif()
