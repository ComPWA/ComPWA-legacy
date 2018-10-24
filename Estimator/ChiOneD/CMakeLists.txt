# Create ChiOneD library.
set(lib_srcs ChiOneD.cpp)
set(lib_headers ChiOneD.hpp)

add_library(ChiOneD
  SHARED ${lib_srcs}  ${lib_headers}
)
target_link_libraries(ChiOneD
  PUBLIC Core
)

#target_include_directories(
#  ChiOneD
#  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
#)

#
# Install
#
install(FILES ${lib_headers}
  DESTINATION include/Estimator/ChiOneD
)
install(TARGETS ChiOneD
  LIBRARY DESTINATION lib
)
