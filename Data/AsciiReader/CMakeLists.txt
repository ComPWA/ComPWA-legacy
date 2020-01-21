# Create AsciiReader library.
set(lib_srcs AsciiReader.cpp)
set(lib_headers AsciiReader.hpp)

add_library( AsciiReader
  ${lib_srcs} ${lib_headers}
)
target_link_libraries(AsciiReader
  Data
)

add_subdirectory(test)

# Install
install(FILES ${lib_headers}
  DESTINATION include/Data/AsciiReader
)
install(TARGETS AsciiReader
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
