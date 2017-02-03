# CMakeLists files in this project can
# refer to the root source directory of the project as ${COMPWA_SOURCE_DIR} and
# to the root binary directory of the project as ${COMPWA_BINARY_DIR}.
cmake_minimum_required (VERSION 3.0)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
project (COMPWA)

#set the default path for built executables to the "bin" directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
#set the default path for built libraries to the "lib" directory
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

set (CMAKE_CXX_STANDARD 11)

#Add Boost
#set(Boost_INCLUDE_DIR $(BOOST_ROOT)/build/include)
#set(Boost_LIBRARY_DIR $(BOOST_LIBRARYDIR))
find_package(Boost COMPONENTS system serialization thread date_time log log_setup REQUIRED)
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
#add_subdirectory (Physics)
#add_subdirectory (executables)
