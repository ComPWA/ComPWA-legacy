# Create AsciiReader library.
set(lib_srcs AsciiReader.cpp)
set(lib_headers AsciiReader.hpp)

add_library( AsciiReader
  SHARED ${lib_srcs} ${lib_headers}
)
target_link_libraries(AsciiReader
  Data
)

#
# Install
#
install(FILES ${lib_headers}
  DESTINATION include/Data/AsciiReader
)
install(TARGETS AsciiReader
  LIBRARY DESTINATION lib
)
