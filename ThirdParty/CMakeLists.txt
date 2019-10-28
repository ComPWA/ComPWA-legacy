# Get system-wide TBB
# System-wide search must be commented out here since pstl's 
# FindTBB.cmake is buggy and does not work if find_package( TBB ) 
# were called before 
#find_package( TBB 2018 )

# If no system-wide version is found, use the local one
if(NOT TBB_FOUND)
  message(STATUS "System-wide TBB not found. Use local one.")
  set(TBB_BUILD_TESTS OFF CACHE BOOL "")
  set(TBB_BUILD_STATIC OFF CACHE BOOL "")
  set(TBB_BUILD_SHARED ON CACHE BOOL "")
  add_subdirectory(tbb)

  # We do not care about compiler warnings of thirdparty dependencies
  target_compile_options(tbb PRIVATE "-w")
  target_compile_options(tbbmalloc PRIVATE "-w")
  target_compile_options(tbbmalloc_proxy PRIVATE "-w")

  # The fork github.com/wjakob/tbb does unfortunatly not properly export 
  # its targets
  add_library(TBB::tbb ALIAS tbb)
  add_library(TBB::tbbmalloc ALIAS tbbmalloc)
  add_library(TBB::tbbmalloc_proxy ALIAS tbbmalloc_proxy)
  export(TARGETS tbb tbbmalloc tbbmalloc_proxy
    NAMESPACE TBB:: FILE ${CMAKE_CURRENT_BINARY_DIR}/TBBConfig.cmake)

  # The export command above is executed after(!) find_package. Therefore we
  # have to create a dummy TBBConfig.cmake file to fool find_library. This 
  # will be overwritten later
  file(WRITE ${CMAKE_BINARY_DIR}/ThirdParty/tbb/TBBConfig.cmake "")
  include(CMakePackageConfigHelpers)
  write_basic_package_version_file(
    ${CMAKE_BINARY_DIR}/ThirdParty/tbb/TBBConfigVersion.cmake 
    VERSION 2019.0
    COMPATIBILITY AnyNewerVersion)

  # We need to set some variables in order to trick the CONFIG mode of pstl's
  # FindTBB.cmake
  set(TBB_DIR ${CMAKE_BINARY_DIR}/ThirdParty/tbb CACHE PATH "")
  set(TBB_IMPORTED_TARGETS "TBB::tbb;TBB::tbbmalloc;TBB::tbbmalloc_proxy")
  set(TBB_tbb_FOUND TRUE)
  set(TBB_tbbmalloc_FOUND TRUE)
  set(TBB_tbbmalloc_proxy_FOUND TRUE)
endif()

add_subdirectory(easyloggingpp)
add_subdirectory(qft++)
add_subdirectory(parallelstl)
add_subdirectory(EvtGen)
