# Create ChiOneD library.
set(SOURCE_LIST ChiOneD.cpp)
set(HEADER_LIST ChiOneD.hpp)

add_library(ChiOneD ${SOURCE_LIST} ${HEADER_LIST})
target_link_libraries(ChiOneD PUBLIC Core)

#
# Install
#
install(FILES ${HEADER_LIST} DESTINATION include/Estimator/ChiOneD)
install(
  TARGETS ChiOneD
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
