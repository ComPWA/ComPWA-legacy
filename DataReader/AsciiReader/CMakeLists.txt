# Create AsciiReader library.
SET( lib_srcs AsciiReader.cpp )
SET( lib_headers AsciiReader.hpp )

add_library( AsciiReader
  SHARED ${lib_srcs} ${lib_headers}
)
target_link_libraries( AsciiReader
  DataReader
)

target_include_directories( AsciiReader
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

install (FILES ${lib_headers} DESTINATION include/ComPWA)
