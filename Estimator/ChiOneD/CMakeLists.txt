# Create ChiOneD library.
set(lib_srcs ChiOneD.cpp)
set(lib_headers ChiOneD.hpp)

add_library(ChiOneD ${lib_srcs} ${lib_headers})
target_link_libraries(ChiOneD PUBLIC Core)

#
# Install
#
install(FILES ${lib_headers} DESTINATION include/Estimator/ChiOneD)
install(
  TARGETS ChiOneD
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
