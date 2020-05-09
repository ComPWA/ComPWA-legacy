# Create MinLogLH library.

set(SOURCE_LIST MinLogLH.cpp SumMinLogLH.cpp)
set(HEADER_LIST MinLogLH.hpp SumMinLogLH.hpp)

add_library(MinLogLH ${SOURCE_LIST} ${HEADER_LIST})

target_link_libraries(
  MinLogLH
  PUBLIC Core
         FunctionTree
         Data
)

install(FILES ${HEADER_LIST}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Estimator/MinLogLH
)
install(
  TARGETS MinLogLH
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
