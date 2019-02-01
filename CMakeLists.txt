#
# ComPWA - A Common framework for Partial Wave Analysis
#
project( COMPWA )

#
# CMake configuration
#
cmake_minimum_required (VERSION 3.4)
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Setting policy CMP0060 to the new version enforces that libraries are linked
# using their full path. That should help in case that multiple versions of a
# library (e.g. boost) are installed on a system
cmake_policy(SET CMP0060 NEW)

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)
# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)



# Which target should be build by default?
set(DEFAULT_BUILD_TYPE "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

message( STATUS "Compiler: " ${CMAKE_CXX_COMPILER_ID} " ( " ${CMAKE_CXX_COMPILER} " )" )
message( STATUS "All possible build types and flags are listed below:")
message( STATUS "DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")
message( STATUS "RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")
message( STATUS "RELWITHDEBINFO is ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
message( STATUS "MINSIZEREL is ${CMAKE_CXX_FLAGS_MINSIZEREL}")
message( STATUS "Use this CMake flag to select build type: -DCMAKE_BUILD_TYPE=BUILD_TYPE")
message( STATUS "Using build type: " ${CMAKE_BUILD_TYPE} )

# Additonal compile flags for various compilers
# ${CMAKE_CXX_COMPILER_ID} can be one of {GNU Clang AppleClang Intel MSVC}
# - verbose output on loop vectorization
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang") # using Clang or AppleClang
    #set(CMAKE_CXX_FLAGS
    #    "${CMAKE_CXX_FLAGS} -Rpass-analysis=loop-vectorize -Rpass=loop-vectorize")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") # using GCC
    # gcc 4.8/4.9 have for example no full regex implementation 
    # and will cause runtime errors
    # full c++11 support is only given in gcc 5.1
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "5.1")
        message(FATAL_ERROR "Version of gcc is too low, and does not have full c++11 support.")
        message(FATAL_ERROR "Please install gcc 5.1 or higher")
    endif()
    #set(CMAKE_CXX_FLAGS
    #    "${CMAKE_CXX_FLAGS} -ftree-vectorize -ftree-vectorizer-verbose=1")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel") # using Intel C++
     # nothing here
endif()

# - Enable all warnings
# This is not good code, in the sense that its not portable (across compilers)
# The portable way using target_compile_options is just so inconvienient, that
# we use this for now. 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

# Enable -fPIC flag
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

message( STATUS "Using CXX compiler flags: " ${CMAKE_CXX_FLAGS} )

# Setting name prefix for libraries
set( CMAKE_SHARED_LIBRARY_PREFIX "libComPWA_" )

#
# External dependencies
#

# Frameworks for parallelization
find_package( OpenMP QUIET )
if (${OPENMP_CXX_FOUND})
  message( STATUS "OpenMP version " ${OpenMP_CXX_VERSION} )
endif()

#
# Boost unit_test_framework and serialization are used throughout the software
# so that they are a requirement
#
find_package(Boost COMPONENTS
	unit_test_framework
    serialization
    filesystem
    system
   	program_options
    REQUIRED
)
# Testing purposes
# set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH}:${Boost_LIBRARY_DIRS}")
# set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH}:/cluster/easybuild/broadwell/software/compiler/GCCcore/6.3.0/lib64")

# Use dynamic libraries for log, test and serialization
#add_definitions(-DBOOST_LOG_DYN_LINK)
#add_definitions(-DBOOST_TEST_DYN_LINK)
#add_definitions(-DBOOST_SERIALIZATION_DYN_LINK)

# find python3
# unfortunately there is a change in the interface of CMake at version 12
# in which the PythonInterp gets deprecated and FindPython3 should be used
set(Python3_FOUND FALSE)
if (${CMAKE_MAJOR_VERSION} EQUAL 3 AND ${CMAKE_MINOR_VERSION} LESS 12)
    find_package(PythonInterp)
    if(${PYTHONINTERP_FOUND} AND ${PYTHON_VERSION_MAJOR} GREATER 2)
        set(Python3_FOUND TRUE)
    endif()
else()
    find_package(Python3)
endif()

find_package(ROOT COMPONENTS Minuit2)

find_package(GSL COMPONENTS gsl gslcblas)

# Geneva minimizer module. Check first if library is available.
find_package(Geneva)

#
# Third party libraries included in the repository
# - easyloggingpp
# - pybind11
# - qft++
#
add_subdirectory(ThirdParty)

# setting some global include directories (this is not the best of practices)
include_directories(${COMPWA_SOURCE_DIR})

#
# Enable target 'test'
#
enable_testing()
set(CTEST_OUTPUT_ON_FAILURE TRUE)
 
#
# Submodules
#
add_subdirectory(Core)
add_subdirectory(Tools)
add_subdirectory(Data)
add_subdirectory(Estimator)
add_subdirectory(Optimizer)
add_subdirectory(Physics)
add_subdirectory(Examples)

#
# Doxygen
#
# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile
                 ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)

    # Adding custom command to generate documenation
    # The file doxygen.stamp is generated to keep track
    # if documentation was generated. Unfortunatly cmake copy 
    # command does not support wildcartds, so we have to copy 
    # file-by-file.
    add_custom_command(
        OUTPUT
            doc/doxygen.stamp
        DEPENDS
            doc/Doxyfile
        COMMAND 
            ${CMAKE_COMMAND} -E make_directory
            ${CMAKE_CURRENT_BINARY_DIR}/doc/images/
        COMMAND 
            ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_SOURCE_DIR}/doc/images/ComPWA_Modules.jpg
            ${CMAKE_CURRENT_BINARY_DIR}/doc/images
        COMMAND 
            ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_SOURCE_DIR}/doc/images/HelicityAngle.png
            ${CMAKE_CURRENT_BINARY_DIR}/doc/images
        COMMAND 
            ${CMAKE_COMMAND} -E make_directory
            ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen-theme/
        COMMAND 
            ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen-theme/header.html
            ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen-theme/
        COMMAND 
            ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen-theme/footer.html
            ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen-theme/
        COMMAND 
            ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen-theme/customdoxygen.css
            ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen-theme/
        COMMAND 
            ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen-theme/DoxygenLayout.xml
            ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen-theme/
        COMMAND 
            ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen-theme/doxy-boot.js
            ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen-theme/
        COMMAND
            doxygen doc/Doxyfile
        COMMAND
            cmake -E touch doc/doxygen.stamp

        COMMENT
            "Generating API documentation with Doxygen"
        VERBATIM
    )

    # Actual doc target depends on doxygen.stamp and runs 
    # custom command if needed
    add_custom_target(
        #doc ALL # run the target with make
        doc # run it separatly
        DEPENDS doc/doxygen.stamp
    )
endif(DOXYGEN_FOUND)
