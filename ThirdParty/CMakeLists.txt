# Get system-wide TBB
find_package( TBB 2018 )
# If no system-wide version is found, use the local one
if(NOT TBB_FOUND)
  message(STATUS "System-wide TBB not found. Link local one statically.")
  # set(TBB_BUILD_STATIC ON)
  set(TBB_BUILD_SHARED OFF) # Link statically
  add_subdirectory(tbb)
  set(TBB_FOUND ON)
  add_library(TBB::tbb ALIAS tbb_static)
  add_library(TBB::tbbmalloc ALIAS tbbmalloc_static)
  add_library(TBB::tbbmalloc_proxy ALIAS tbbmalloc_proxy_static)
  export(TARGETS tbb_static tbbmalloc_static tbbmalloc_proxy_static
    NAMESPACE TBB:: FILE TBBTargets.cmake)
  export(PACKAGE TBB)
  # Update version if submodule is updates
  set(TBB_VERSION "2019(U2)") 
  set(TBB_INTERFACE_VERSION "11002")
endif()

add_subdirectory(easyloggingpp)
add_subdirectory(pybind11)
add_subdirectory(qft++)
add_subdirectory(parallelstl)
add_subdirectory(EvtGen)
