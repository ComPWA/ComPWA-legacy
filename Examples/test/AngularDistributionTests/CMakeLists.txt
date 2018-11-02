if( ${Python3_FOUND}
    AND TARGET Minuit2IF AND TARGET Data AND TARGET HelicityFormalism
    AND TARGET Tools AND TARGET pycompwa
)

set(TEST_NAMES
  D1ToD0PipPim
)

foreach(TestName ${TEST_NAMES})
  set(CTestName "Examples_AngularDistributionTests_${TestName}")

  add_test(NAME ${CTestName}
    COMMAND pytest -v -k ${TestName}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
  set_property(TEST ${CTestName} PROPERTY
    ENVIRONMENT "PYTHONPATH=${COMPWA_SOURCE_DIR}/Physics/ExpertSystem:
      ${COMPWA_SOURCE_DIR}/Tools:${CMAKE_BINARY_DIR}/Tools/PythonInterface")

endforeach()

else ()
  message(WARNING "Requirements not found! Not building Examples AngularDistribution Tests!")
endif()
