# Require external version of tbb
# find_package( TBB 2018 REQUIRED )

# The problem with the internal version is that no install target 
# is provided and currently I don't know how to work around that. Furthermore,
# there are incompatibilities between the FindTBB.cmake and the proposed 
# installation procedure of tbb.
# if(NOT ${TBB_FOUND})
  # message( STATUS "TBB not found. Build the one included in ComPWA...")
  set(TBB_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/tbb)
  include(${TBB_ROOT}/cmake/TBBBuild.cmake)
  tbb_build(TBB_ROOT ${TBB_ROOT} CONFIG_DIR TBB_DIR)
  # tbb_build(TBB_ROOT ${TBB_ROOT}
            # MAKE_ARGS "tbb_build_prefix=tbb;tbb_build_dir=${CMAKE_CURRENT_BINARY_DIR}/tbb")
# endif()

add_subdirectory(easyloggingpp)
#add_subdirectory(backward-cpp)
add_subdirectory(pybind11)
add_subdirectory(qft++)
add_subdirectory(parallelstl)
message(STATUS "Parallel STL uses TBB at ${TBB_LIBRARIES}")
add_subdirectory(EvtGen)
