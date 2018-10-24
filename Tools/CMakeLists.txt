#################
# Tools Library #
#################

# Check if all requirements are found. Skip modules otherwise.
if( ${ROOT_FOUND} AND ${GSL_FOUND}
#    AND TARGET DataReader AND TARGET HelicityFormalism
)

FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)

### Generator dictionaries
# Currently we do not generate dictionaries of compatibility reasons
# We have to add this directory to the include directories
# so rootcint can find the header and LinkDef files
#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
#ROOT_GENERATE_DICTIONARY( G__Dict SimpleEfficiency.hpp LINKDEF libDPKinematicsLinkDef.h )
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/G__Dict.cxx.o DESTINATION lib/)
#add_library ( DPKinematics SHARED ${lib_srcs} G__Dict.cxx ${lib_headers} )

add_library(Tools
  SHARED ${lib_srcs} ${lib_headers}
)

target_link_libraries(Tools
  PUBLIC
  	Core
  	DataReader
  	HelicityFormalism
  	ROOT::Physics
  	ROOT::MathCore
  	ROOT::Core
  PRIVATE
  	ROOT::EG
  	ROOT::Tree
  	ROOT::Hist
  	ROOT::Gpad
  	ROOT::Graf
  	ROOT::RIO
  INTERFACE
  	${GSL_LIBRARIES}
)

target_include_directories(Tools
  PUBLIC 
    $<INSTALL_INTERFACE:include>
    ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} ${GSL_INCLUDE_DIR}
)

#
# Install
#
install(FILES ${lib_headers}
    DESTINATION include/Tools
)
install(TARGETS Tools
    LIBRARY DESTINATION lib
)

#
# TESTING
#
# Testing routines are stored in separate directory

add_executable(toolsTest test/ToolsTest.cpp)

target_link_libraries(toolsTest
    Core
    Boost::unit_test_framework
)

set_target_properties(toolsTest
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
)
target_compile_definitions(toolsTest PRIVATE "BOOST_TEST_DYN_LINK=1")

add_test(NAME toolsTest
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/toolsTest
)

else ()
  message(WARNING "Requirements not found! Not building Tools module!")
endif()


add_subdirectory(PythonInterface)
add_subdirectory(Plotting)
