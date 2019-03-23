# Get system-wide TBB
# find_package( TBB 2018 )
# If no system-wide version is found, use the local one
if(NOT TBB_FOUND)
  message(STATUS "System-wide TBB not found. Link local one statically.")
	set(TBB_BUILD_STATIC OFF)
	set(TBB_BUILD_TESTS OFF)
  # set(TBB_BUILD_SHARED OFF) # Link statically
	add_subdirectory(tbb)
  # include(ExternalProject)
  # ExternalProject_Add (
  #     TBB
  #     # GIT_REPOSITORY "https://github.com/wjakob/tbb.git"
  #     # GIT_TAG "b066def"
	#     SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tbb"
  #     UPDATE_COMMAND ""
  #     INSTALL_COMMAND "" 
  #     CMAKE_CACHE_ARGS "-DTBB_BUILD_STATIC:BOOL=OFF -DTBB_BUILD_TESTS:BOOL=OFF"
	# )
	# set(TBB_DIR ${CMAKE_BINARY_DIR}/ThirdParty/tbb)
	set(TBB_INSTALL_DIR ${CMAKE_BINARY_DIR}/ThirdParty/tbb)
	set(TBB_DIR ${CMAKE_BINARY_DIR}/ThirdParty/tbb)
  find_package(TBB CONFIG )
  message(STATUS "---> ${TBB_FOUND} - ${TBB_DIR}")

  # set(TBB_FOUND ON)
  # set(TBB_DIR ${CURRENT_SOURCE_DIR}/tbb)
	# add_library(TBB::tbb ALIAS tbb)
	# add_library(TBB::tbbmalloc ALIAS tbbmalloc)
	# add_library(TBB::tbbmalloc_proxy ALIAS tbbmalloc_proxy)
	# export(TARGETS tbb_static tbbmalloc_static tbbmalloc_proxy_static
		# FILE ${PROJECT_BINARY_DIR}/ThirdParty/TBBConfig.cmake)
  # Update version if submodule is updates
  # set(TBB_VERSION "2019(U2)")
  # set(TBB_INTERFACE_VERSION "11002")

  # ... for the build tree
  # set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
  # configure_file(TBBConfig.cmake.in
  #   "${PROJECT_BINARY_DIR}/ThirdParty/TBBConfig.cmake" @ONLY)
  # # ... for the install tree
  # set(CONF_INCLUDE_DIRS "\${TBB_CMAKE_DIR}/include")
  # configure_file(TBBConfig.cmake.in
  #   "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/TBBConfig.cmake" @ONLY)
  # # ... for both
  # configure_file(TBBConfigVersion.cmake.in
  #   "${PROJECT_BINARY_DIR}/TBBConfigVersion.cmake" @ONLY)
  #
  # # Install the TBBConfig.cmake and TBBConfigVersion.cmake
  # install(FILES
  #   "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/TBBConfig.cmake"
  #   "${PROJECT_BINARY_DIR}/TBBConfigVersion.cmake"
  #   DESTINATION "lib/cmake/TBB" COMPONENT dev)
  #
  # # Install the export set for use with the install-tree
  # install(EXPORT TBBTargets DESTINATION
  #   "lib/cmake/TBB" COMPONENT dev)

endif()

add_subdirectory(easyloggingpp)
add_subdirectory(pybind11)
add_subdirectory(qft++)
# add_subdirectory(parallelstl)
add_subdirectory(EvtGen)
