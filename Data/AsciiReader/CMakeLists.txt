# Create AsciiReader library.
set(lib_srcs AsciiReader.cpp)
set(lib_headers AsciiReader.hpp)

add_library( AsciiReader
  SHARED ${lib_srcs} ${lib_headers}
)
target_link_libraries(AsciiReader
  Data
)
target_include_directories(AsciiReader
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
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
