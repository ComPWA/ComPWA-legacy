# Create MinLogLH library.

set(lib_srcs MinLogLH.cpp SumMinLogLH.cpp)
set(lib_headers MinLogLH.hpp SumMinLogLH.hpp)

add_library(MinLogLH ${lib_srcs} ${lib_headers})

target_link_libraries(MinLogLH PUBLIC Core FunctionTree Data)

install(FILES ${lib_headers} DESTINATION include/Estimator/MinLogLH)
install(
  TARGETS MinLogLH
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
