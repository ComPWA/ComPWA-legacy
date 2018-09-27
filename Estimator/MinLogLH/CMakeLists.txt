# Create MinLogLH library.

set(lib_srcs MinLogLH.cpp SumMinLogLH.cpp)
set(lib_headers MinLogLH.hpp SumMinLogLH.hpp)

add_library(MinLogLH
  SHARED ${lib_srcs}
)

target_link_libraries(MinLogLH
  PUBLIC Core
)

#
# Install
#
install(FILES ${lib_headers}
  DESTINATION include/Estimator/MinLogLH
)
install(TARGETS MinLogLH
  LIBRARY DESTINATION lib
)