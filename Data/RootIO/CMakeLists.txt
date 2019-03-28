if(${ROOT_FOUND})

set(lib_srcs
  RootDataIO.cpp
  RootEfficiency.cpp
)

set(lib_headers
  RootDataIO.hpp
  RootEfficiency.hpp
)

add_library(RootDataIO
  SHARED ${lib_srcs} ${lib_headers}
)

target_link_libraries(RootDataIO
  PUBLIC Core Data ROOT::Hist ROOT::Core
  PRIVATE ROOT::EG ROOT::Physics ROOT::Tree ROOT::RIO
)

# target_include_directories(RootDataIO
  # PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR}
# )

#
# Install
#
install (FILES ${lib_headers}
  DESTINATION include/Data/RootDataIO
)
install(TARGETS RootDataIO
	LIBRARY DESTINATION lib
)

#
# TESTING
#
# Testing routines are stored in separate directory

add_executable(Data_RootDataIOTest test/RootDataIOTest.cpp)

target_link_libraries(Data_RootDataIOTest
  RootDataIO
  Tools
  Boost::unit_test_framework
)

set_target_properties(Data_RootDataIOTest
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
)

add_test(NAME Data_RootDataIOTest
	WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/Data_RootDataIOTest
)

else ()
	message(WARNING "ROOT not found! Not building RootDataIO!")
endif()
