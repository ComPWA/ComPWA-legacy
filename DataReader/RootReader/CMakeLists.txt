#############################
# Create RootReader library #
#############################

if(${ROOT_FOUND})

set(lib_srcs
  RootReader.cpp
  RootEfficiency.cpp
)

set(lib_headers
  RootReader.hpp
  RootEfficiency.hpp
)

add_library(RootReader
  SHARED ${lib_srcs} ${lib_headers}
)

target_link_libraries(RootReader
  PUBLIC Core DataReader ROOT::Hist ROOT::Core
  PRIVATE ROOT::EG ROOT::Physics ROOT::Tree ROOT::RIO
)

target_include_directories(RootReader
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR}
)

#
# Install
#
install (FILES ${lib_headers}
  DESTINATION include/DataReader/RootReader
)
install(TARGETS RootReader
	LIBRARY DESTINATION lib
)

#
# TESTING
#
# Testing routines are stored in separate directory

add_executable(RootReader_RootReaderTest test/RootReaderTest.cpp)

target_link_libraries(RootReader_RootReaderTest
  RootReader
  Tools
  Boost::unit_test_framework
)

set_target_properties(RootReader_RootReaderTest
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
)

target_compile_definitions(RootReader_RootReaderTest PRIVATE "BOOST_TEST_DYN_LINK=1")

add_test(NAME RootReader_RootReaderTest
	WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/RootReader_RootReaderTest
)

else ()
	message(WARNING "ROOT not found! Not building RootReader!")
endif()
