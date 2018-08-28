#
# ComPWA - A Common framework for Partial Wave Analysis
#
project( COMPWA )

#
# CMake configuration
#
cmake_minimum_required (VERSION 3.4)
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Setting policy CMP0060 to the new version enforces that 
# libraries are linked using their full path. That should help in
# case that multiple versions of a library (e.g. boost) are
# installed on a system
cmake_policy(SET CMP0060 NEW)

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

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)

message( STATUS "Compiler: " ${CMAKE_CXX_COMPILER_ID} " ( " ${CMAKE_CXX_COMPILER} " )" )
message( STATUS "All possible build types and flags are listed below:")
message( STATUS "DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")
message( STATUS "RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")
message( STATUS "RELWITHDEBINFO is ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
message( STATUS "MINSIZEREL is ${CMAKE_CXX_FLAGS_MINSIZEREL}")
message( STATUS "Use this CMake flag to select build type: -DCMAKE_BUILD_TYPE=BUILD_TYPE")
message( STATUS "Using build type: " ${CMAKE_BUILD_TYPE} )

#
# Additonal compile flags for various compilers
# ${CMAKE_CXX_COMPILER_ID} can be one of {GNU Clang AppleClang Intel MSVC}
# - verbose output on loop vectorization
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang") # using Clang or AppleClang
    #set(CMAKE_CXX_FLAGS
    #    "${CMAKE_CXX_FLAGS} -Rpass-analysis=loop-vectorize -Rpass=loop-vectorize")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") # using GCC
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
find_package( TBB )
if (${TBB_FOUND})
  message( STATUS "TBB found version ${TBB_VERSION}" )
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
# Use dynamic libraries for log, test and serialization
add_definitions(-DBOOST_LOG_DYN_LINK)
add_definitions(-DBOOST_TEST_DYN_LINK)
add_definitions(-DBOOST_SERIALIZATION_DYN_LINK)

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
include_directories(${ELPP_INCLUDE_DIRS})


#
# Enable target 'test'
#
enable_testing()
set(CTEST_OUTPUT_ON_FAILURE TRUE)

# also run pytest for testing
add_test(NAME PythonTests
         COMMAND pytest -v 
         --ignore=Tools/PythonInterface
         --ignore=ThirdParty
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
set_property(TEST PythonTests PROPERTY
    ENVIRONMENT "PYTHONPATH=${CMAKE_SOURCE_DIR}/Physics/ExpertSystem")
  
#
# Submodules
#
add_subdirectory(Core)
add_subdirectory(DataReader)
add_subdirectory(Estimator)
add_subdirectory(Optimizer)
add_subdirectory(Physics)
add_subdirectory(Tools)
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

    # Old solution for doxygen
    #add_custom_target(doc
          #${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
          #WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
          #COMMENT "Generating API documentation with Doxygen" VERBATIM)

    # Copy images
    #file(GLOB img
    #    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/doc/images *
    #)

    #add_custom_command(
      #TARGET doc PRE_BUILD
      #COMMAND ${CMAKE_COMMAND} -E make_directory
      #${CMAKE_CURRENT_BINARY_DIR}/doc/images/ )
    #add_custom_command(
      #TARGET doc PRE_BUILD
      #COMMAND ${CMAKE_COMMAND} -E make_directory
      #${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen-theme/ )
    #add_custom_command(
      #TARGET doc PRE_BUILD
      #COMMAND ${CMAKE_COMMAND} -E copy
      #${CMAKE_CURRENT_SOURCE_DIR}/doc/images/*
      #${CMAKE_CURRENT_BINARY_DIR}/doc/images/ )
    #add_custom_command(
      #TARGET doc PRE_BUILD
      #COMMAND ${CMAKE_COMMAND} -E copy
      #${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen-theme/*
      #${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen-theme/ )

endif(DOXYGEN_FOUND)
