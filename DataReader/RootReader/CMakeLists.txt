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
  SHARED ${lib_srcs}
)

target_link_libraries(RootReader
  PUBLIC Core DataReader ROOT::Core ROOT::Hist
  PRIVATE ROOT::EG ROOT::Physics ROOT::Tree ROOT::RIO
)

target_include_directories(RootReader
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR}
)

#
# Install
#
install (FILES ${lib_headers}
    DESTINATION include/ComPWA/DataReader/RootReader
)

install(TARGETS RootReader
	LIBRARY DESTINATION lib/ComPWA
)

#
# TESTING
#
# Testing routines are stored in separate directory

add_executable(rootReaderTest test/RootReaderTest.cpp)

target_link_libraries(rootReaderTest
    RootReader
    Tools
    Boost::unit_test_framework
)

set_target_properties(rootReaderTest
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
)

add_test(NAME RootReaderTest
	WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/rootReaderTest
)

else ()
	message(WARNING "ROOT not found! Not building RootReader!")
endif()
