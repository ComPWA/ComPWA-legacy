#############################
# Create DataReader library #
#############################

set(lib_srcs
  Data.cpp
  DataCorrection.cpp
  CorrectionTable.cpp
)
set(lib_headers
  Data.hpp
  DataCorrection.hpp
  CorrectionTable.hpp
)

add_library(DataReader
  SHARED ${lib_srcs}
)

target_link_libraries(DataReader
  PUBLIC Core
)

target_include_directories(DataReader
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

add_subdirectory(AsciiReader)
add_subdirectory(RootReader)

#
# Install
#
install (FILES ${lib_headers}
    DESTINATION include/ComPWA/DataReader
)
install(TARGETS DataReader
	LIBRARY DESTINATION lib/ComPWA
)
