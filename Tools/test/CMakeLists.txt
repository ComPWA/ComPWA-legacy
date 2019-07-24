# Check if all requirements are found
if(TARGET Tools)

add_executable(UpdatePTreeParameterTest UpdatePTreeParameterTest.cpp)
target_link_libraries(UpdatePTreeParameterTest
    Tools
    Boost::unit_test_framework
)
set_target_properties(UpdatePTreeParameterTest
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
)

add_test(NAME UpdatePTreeParameterTest
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/UpdatePTreeParameterTest
)

add_executable(SaveAndLoadFitResultTest SaveAndLoadFitResultTest.cpp)

target_link_libraries(SaveAndLoadFitResultTest
    Core RootData Minuit2IF MinLogLH Tools HelicityFormalism
    ROOT::Core
    Boost::unit_test_framework
)
set_target_properties(SaveAndLoadFitResultTest
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
)

add_test(NAME SaveAndLoadFitResultTest
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/SaveAndLoadFitResultTest
)
endif()
