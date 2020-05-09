# Project definition
cmake_minimum_required(VERSION 3.4 FATAL_ERROR)
project(
  COMPWA
  DESCRIPTION "The Common Partial Wave Analysis framework"
  LANGUAGES CXX
)
set(CMAKE_CXX_STANDARD 14)

# CMake configuration
list(
  APPEND
  CMAKE_MODULE_PATH
  "${CMAKE_SOURCE_DIR}/cmake/Modules/"
)

# Which target should be build by default?
set(DEFAULT_BUILD_TYPE "Release")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# CMake policies
# https://cmake.org/cmake/help/latest/policy/CMP0060.html
cmake_policy(
  SET
  CMP0060
  NEW
)
if(POLICY CMP0048)
  # https://cmake.org/cmake/help/latest/policy/CMP0048.html
  cmake_policy(
    SET
    CMP0048
    NEW
  )
endif()

# Configure RPATH
if(APPLE)
  list(
    APPEND
    CMAKE_INSTALL_RPATH
    "@loader_path/../lib;@loader_path"
  )
else()
  list(
    APPEND
    CMAKE_INSTALL_RPATH
    "$ORIGIN/:$ORIGIN/../lib"
  )
endif()
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(BUILD_SHARED_LIBS ON)

# Setting name prefix for libraries
set(CMAKE_SHARED_LIBRARY_PREFIX "libComPWA_")

# Build options
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(
    STATUS
      "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified."
  )
  set(CMAKE_BUILD_TYPE
      "${DEFAULT_BUILD_TYPE}"
      CACHE STRING
            "Choose the type of build."
            FORCE
  )
  # Set the possible values of build type for cmake-gui
  set_property(
    CACHE CMAKE_BUILD_TYPE
    PROPERTY STRINGS
             "Debug"
             "Release"
             "MinSizeRel"
             "RelWithDebInfo"
  )
endif()

# Compiler flags
if("${CMAKE_CXX_COMPILER_ID}"
   STREQUAL
   "GNU"
) # using GCC
  if(CMAKE_CXX_COMPILER_VERSION
     VERSION_LESS
     "5.1"
  )
    message(FATAL_ERROR "Version of gcc is too low, and does not have \
    full c++11 support. Please install gcc 5.1 or higher."
    )
  endif()
endif()
message(STATUS "Global CXX compiler flags: " ${CMAKE_CXX_FLAGS})

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Code Coverage Configuration
option(
  CODE_COVERAGE
  "Enable coverage reporting"
  OFF
)
if(CODE_COVERAGE
   AND CMAKE_CXX_COMPILER_ID
       MATCHES
       "GNU|Clang"
)
  set(CMAKE_CXX_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
endif()

option(
  SAFE_LINKING
  "Enable to avoid conflicts between different versions \
of libraries stdc++, boost, ...)"
  OFF
)
if(SAFE_LINKING)
  if(NOT APPLE)
    # Setting RPATH instead of RUNPATH
    # RPATH is search before LD_LIBRARY_PATH at runtime
    set(CMAKE_EXE_LINKER_FLAGS
        "${CMAKE_EXE_LINKER_FLAGS} -Wl,--disable-new-dtags"
    )
    set(CMAKE_SHARED_LINKER_FLAGS
        "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--disable-new-dtags"
    )
  endif()
  if(NOT
     "${CMAKE_CXX_COMPILER_ID}"
     MATCHES
     "Clang"
  )
    message(STATUS "Option SAFE_LINKING includes static linking of stdc++ "
                   "GCC."
    )
    set(CMAKE_EXE_LINKER_FLAGS
        "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++ -static-libgcc"
    )
    set(CMAKE_SHARED_LINKER_FLAGS
        "${CMAKE_SHARED_LINKER_FLAGS} -static-libstdc++ -static-libgcc"
    )
  endif()
endif()

# External dependencies
if(SAFE_LINKING)
  message(STATUS "Option SAFE_LINKING is enabled. Make sure static boost "
                 "libraries are available and are compiled with -fPIC."
  )
  set(Boost_USE_STATIC_LIBS ON)
else()
  add_definitions(-DBOOST_TEST_DYN_LINK=1)
endif()
find_package(
  Boost
  COMPONENTS unit_test_framework
             serialization
             program_options
             filesystem
  REQUIRED
)
set(THREADS_PREFER_PTHREAD_FLAG ON)

# Add external libraries
find_package(Threads REQUIRED)
find_package(ROOT COMPONENTS Minuit2)
option(
  USE_GENEVA
  "Switch Geneva On/Off"
  OFF
)
if(USE_GENEVA)
  find_package(Geneva QUIET)
endif()
add_subdirectory(ThirdParty)

# Setting ComPWA source dir as include directory
include_directories(${COMPWA_SOURCE_DIR})

# Enable target 'test'
enable_testing()
set(CTEST_OUTPUT_ON_FAILURE TRUE)

# Submodules
add_subdirectory(Core)
add_subdirectory(Tools)
add_subdirectory(Data)
add_subdirectory(Estimator)
add_subdirectory(Optimizer)
add_subdirectory(Physics)
add_subdirectory(Examples)
