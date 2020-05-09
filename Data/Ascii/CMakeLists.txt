set(SOURCE_LIST AsciiDataIO.cpp AsciiHeaderIO.cpp)
set(HEADER_LIST AsciiDataIO.hpp AsciiHeaderIO.hpp)
add_library(Ascii ${SOURCE_LIST} ${HEADER_LIST})
target_link_libraries(Ascii Data)

add_subdirectory(test)

# Install
install(FILES ${HEADER_LIST} DESTINATION include/Data/Ascii)
install(
  TARGETS Ascii
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
