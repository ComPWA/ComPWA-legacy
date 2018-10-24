#############################
#    Create Data library    #
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
  SHARED ${lib_srcs} ${lib_headers}
)

target_link_libraries(Data
  PUBLIC Core
)

target_include_directories(Data
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

add_subdirectory(AsciiReader)
add_subdirectory(RootReader)

#
# Install
#
install (FILES ${lib_headers}
    DESTINATION include/Data
)
install(TARGETS Data
	LIBRARY DESTINATION lib
)
