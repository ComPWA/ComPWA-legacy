# Create AsciiReader library.
SET( lib_srcs AsciiReader.cpp )
SET( lib_headers AsciiReader.hpp )

add_library( ComPWA_AsciiReader
  SHARED ${lib_srcs} ${lib_headers}
)
target_link_libraries( ComPWA_AsciiReader
  ComPWA_DataReader
)

target_include_directories( ComPWA_AsciiReader
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

install (FILES ${lib_headers} DESTINATION include/ComPWA)
