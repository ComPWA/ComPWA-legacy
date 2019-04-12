
# Check if all requirements are found
if( ${ROOT_FOUND} AND TARGET Tools
)

add_executable(RootGeneratorTest RootGeneratorTest.cpp)
target_link_libraries(RootGeneratorTest
    Tools Boost::unit_test_framework
)
set_target_properties(RootGeneratorTest
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
)
add_test(NAME RootGeneratorTest
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/RootGeneratorTest
)
endif()

# Check if all requirements are found
if(TARGET Tools)

add_executable(EvtGenGeneratorTest EvtGenGeneratorTest.cpp)
target_link_libraries(EvtGenGeneratorTest
    Tools
    Boost::unit_test_framework
)
set_target_properties(EvtGenGeneratorTest
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
)

add_test(NAME EvtGenGeneratorTest
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/EvtGenGeneratorTest
)
endif()

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
endif()

# Check if all requirements are found
if (TARGET Tools)

add_executable(SaveAndLoadFitResultTest SaveAndLoadFitResultTest.cpp)

target_link_libraries(SaveAndLoadFitResultTest
    Core Minuit2IF MinLogLH Tools HelicityFormalism
    ${ROOT_LIBARARY} ${Boost_LIBRARY}
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
