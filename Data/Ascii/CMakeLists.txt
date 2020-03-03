set(lib_srcs AsciiDataIO.cpp AsciiHeaderIO.cpp)
set(lib_headers AsciiDataIO.hpp AsciiHeaderIO.hpp)
add_library(Ascii ${lib_srcs} ${lib_headers})
target_link_libraries(Ascii Data)

add_subdirectory(test)

# Install
install(FILES ${lib_headers} DESTINATION include/Data/Ascii)
install(
  TARGETS Ascii
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
