set(lib_srcs
  DataSet.cpp
  DataCorrection.cpp
  CorrectionTable.cpp
)
set(lib_headers
  DataSet.hpp
  DataCorrection.hpp
  CorrectionTable.hpp
)

add_library(Data
  SHARED ${lib_srcs} ${lib_headers}
)

target_link_libraries(Data
  PUBLIC Core
)

target_include_directories(Data
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

install (FILES ${lib_headers}
    DESTINATION include/Data
)
install(TARGETS Data
	LIBRARY DESTINATION lib
)

add_subdirectory(AsciiReader)
add_subdirectory(RootIO)
