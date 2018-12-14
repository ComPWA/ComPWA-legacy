if(TARGET Tools)
#    add_executable(AmpModelTest AmpModel.cpp)
#
#    target_link_libraries(AmpModelTest
#      Core
#      HelicityFormalism
#      Tools
#      Boost::unit_test_framework
#    )
#  
#    target_include_directories(AmpModelTest
#      PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR})
#
#    # Move testing binaries into a testBin directory
#    set_target_properties(AmpModelTest
#      PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
#    )
#    target_compile_definitions(AmpModelTest PRIVATE "BOOST_TEST_DYN_LINK=1")
#
#    # Copy input files for test program
#    add_custom_command(
#        TARGET AmpModelTest POST_BUILD
#        COMMAND ${CMAKE_COMMAND} -E copy
#        ${CMAKE_CURRENT_SOURCE_DIR}/AmpModel-input.xml
#        ${PROJECT_BINARY_DIR}/bin/test/AmpModel-input.xml)
#
#    add_test(NAME AmpModelTest
#      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
#      COMMAND ${PROJECT_BINARY_DIR}/bin/test/AmpModelTest)
      
    # -------------------- Helicity Angle Test -------------------- #      
    add_executable(HelicityAngle HelicityAngle.cpp)

    target_link_libraries(HelicityAngle
      Core
      HelicityFormalism
      Tools
      Boost::unit_test_framework
    )
  
    target_include_directories(HelicityAngle
      PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR})

    # Move testing binaries into a testBin directory
    set_target_properties(HelicityAngle
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    )
    target_compile_definitions(HelicityAngle PRIVATE "BOOST_TEST_DYN_LINK=1")

    add_test(NAME HelicityAngle
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
      COMMAND ${PROJECT_BINARY_DIR}/bin/test/HelicityAngle)

	# -------------------- Helicity Angles Correctness Test -------------------- # 
    add_executable(HelicityAnglesCorrectnessTest HelicityAnglesCorrectnessTest.cpp)

    target_link_libraries(HelicityAnglesCorrectnessTest
      Core
      HelicityFormalism
      Tools
      Boost::unit_test_framework
    )
  
    target_include_directories(HelicityAnglesCorrectnessTest
      PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR})

    # Move testing binaries into a testBin directory
    set_target_properties(HelicityAnglesCorrectnessTest
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    )
    target_compile_definitions(HelicityAnglesCorrectnessTest PRIVATE "BOOST_TEST_DYN_LINK=1")

    # Copy input files for test program
    add_custom_command(
        TARGET HelicityAnglesCorrectnessTest POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_SOURCE_DIR}/HelicityAnglesCorrectnessTest-input.xml
        ${PROJECT_BINARY_DIR}/bin/test/HelicityAnglesCorrectnessTest-input.xml)

    add_test(NAME HelicityAnglesCorrectnessTest
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
      COMMAND ${PROJECT_BINARY_DIR}/bin/test/HelicityAnglesCorrectnessTest)
    
    # -------------------- Helicity Kinematics Tests -------------------- # 
    add_executable(HelicityKinematicsTests HelicityKinematicsTests.cpp)

    target_link_libraries(HelicityKinematicsTests
      Core
      HelicityFormalism
      Boost::unit_test_framework
    )
  
    target_include_directories(HelicityKinematicsTests
      PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR})

    # Move testing binaries into a testBin directory
    set_target_properties(HelicityKinematicsTests
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    )
    target_compile_definitions(HelicityKinematicsTests PRIVATE "BOOST_TEST_DYN_LINK=1")

    add_test(NAME HelicityKinematicsTests
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
      COMMAND ${PROJECT_BINARY_DIR}/bin/test/HelicityKinematicsTests)
else()
  message(WARNING "Requirements not found! Not building tests!")
endif()
