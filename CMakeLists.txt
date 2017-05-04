# CMakeLists files in this project can
# refer to the root source directory of the project as ${COMPWA_SOURCE_DIR} and
# to the root binary directory of the project as ${COMPWA_BINARY_DIR}.
cmake_minimum_required (VERSION 3.0)
#set(CMAKE_BUILD_TYPE Debug) 
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

project( COMPWA )

#set the default path for built executables to the "bin" directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
#set the default path for built libraries to the "lib" directory
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

set (CMAKE_CXX_STANDARD 11)

#using python : 2.7 : /usr/bin : /usr/include/python2.7 : /usr/lib64/python2.7 ;
#using python : 3.5 : /usr/bin/python3 : /usr/include/python3.5m : /usr/lib64/python3.5 ;
#find_package(PythonLibs 2.7 EXACT)
#message(${PYTHON_LIBRARIES})
#INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
#LINK_DIRECTORIES(${PYTHON_LIBRARIES})

#Add Boost
#find_package(Boost COMPONENTS python system filesystem serialization thread chrono date_time log log_setup program_options regex REQUIRED)
find_package(Boost COMPONENTS system filesystem serialization thread chrono date_time log log_setup program_options regex REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

INCLUDE_DIRECTORIES(${COMPWA_SOURCE_DIR})
#LINK_DIRECTORIES(${COMPWA_SOURCE_DIR}/lib)

add_definitions(-DBOOST_LOG_DYN_LINK)
add_definitions(-DBOOST_SERIALIZATION_DYN_LINK)

# Recurse into subdirectories.
add_subdirectory (Core)
add_subdirectory (DataReader)
add_subdirectory (Estimator)
add_subdirectory (Optimizer)
add_subdirectory (Physics)
add_subdirectory (executables)

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)
add_custom_target(doc
${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)
