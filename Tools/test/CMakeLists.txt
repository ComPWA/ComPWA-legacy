
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
