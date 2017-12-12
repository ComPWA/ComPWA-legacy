# CMakeLists files in this project can
# refer to the root source directory of the project as ${COMPWA_SOURCE_DIR} and
# to the root binary directory of the project as ${COMPWA_BINARY_DIR}.
cmake_minimum_required (VERSION 3.0)
#set (CMAKE_BUILD_TYPE Debug)
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")


set (CMAKE_CXX_STANDARD 11)

if (NOT "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  set (CMAKE_CXX_STANDARD 14)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()

# Workaround for clang+ c++14/pybind11
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
 # set (CMAKE_CXX_STANDARD 11)
  set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -std=c++0x")
endif()

#
# ComPWA
#
project( COMPWA )

# set the default path for built executables to the "bin" directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

# set the default path for built libraries to the "lib" directory
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

INCLUDE_DIRECTORIES(${COMPWA_SOURCE_DIR})
#LINK_DIRECTORIES(${COMPWA_SOURCE_DIR}/lib)

# Setting prefix for libraries
set( CMAKE_SHARED_LIBRARY_PREFIX "libComPWA_" )


#
# Frameworks for parallelization
#
find_package( OpenMP QUIET )
if (${OPENMP_CXX_FOUND})
  message( STATUS "OpenMP version " ${OpenMP_CXX_VERSION} )
endif()

find_package( TBB )
if (${TBB_FOUND})
  message( STATUS "TBB found version ${TBB_VERSION}" )
endif()

#
# ROOT installation
#
find_package( ROOT )
find_package( Minuit2 )


#
# FindBoost. unit_test_framework and serialization
# throughout the software, so we require them here
# local dependencies should be required in each
# module
#
find_package( Boost COMPONENTS
               unit_test_framework
               serialization
               REQUIRED
)

#Find Python
#find_package( PythonLibs )

# Workaround for clang+boost.python
# ${CMAKE_CXX_COMPILER_ID} can be one of {GNU Clang AppleClang Intel MSVC}
# The python interface is still work in progress, so we disable python by
# default.
if (NOT "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  #find_package(PythonLibs 2.7 EXACT)
  #INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
  #LINK_DIRECTORIES(${PYTHON_LIBRARIES})
  #find_package( Boost COMPONENTS
  #             python
  #             REQUIRED
  #)
endif()

# Use dynamic libraries for log, test and serialization
add_definitions(-DBOOST_LOG_DYN_LINK)
add_definitions(-DBOOST_TEST_DYN_LINK)
add_definitions(-DBOOST_SERIALIZATION_DYN_LINK)

include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

#
# Enable target 'test'
#
enable_testing()


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
add_subdirectory(executables)

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
