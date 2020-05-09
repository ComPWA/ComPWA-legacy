set(SOURCE_LIST DataSet.cpp Generate.cpp)
set(HEADER_LIST DataSet.hpp Generate.hpp)

add_library(Data ${SOURCE_LIST} ${HEADER_LIST})

target_link_libraries(Data PUBLIC Core pstl::ParallelSTL)

install(
  TARGETS Data
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(FILES ${HEADER_LIST} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Data)

add_subdirectory(Ascii)
add_subdirectory(EvtGen)
add_subdirectory(Root)
