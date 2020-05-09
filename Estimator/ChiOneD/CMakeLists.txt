# Create ChiOneD library.
set(SOURCE_LIST ChiOneD.cpp)
set(HEADER_LIST ChiOneD.hpp)

add_library(ChiOneD ${SOURCE_LIST} ${HEADER_LIST})
target_link_libraries(ChiOneD PUBLIC Core)

#
# Install
#
install(FILES ${HEADER_LIST}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Estimator/ChiOneD
)
install(
  TARGETS ChiOneD
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
