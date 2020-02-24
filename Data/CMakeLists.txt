set(lib_srcs DataSet.cpp Generate.cpp)
set(lib_headers DataSet.hpp Generate.hpp)

add_library(Data ${lib_srcs} ${lib_headers})

target_link_libraries(Data PUBLIC Core pstl::ParallelSTL)

install(
  TARGETS Data
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(FILES ${lib_headers} DESTINATION include/Data)

add_subdirectory(Ascii)
add_subdirectory(EvtGen)
add_subdirectory(Root)
