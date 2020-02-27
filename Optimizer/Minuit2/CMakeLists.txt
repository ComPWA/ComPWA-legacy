set(lib_srcs MinuitIF.cpp MinuitResult.cpp)
set(lib_headers MinuitFcn.hpp MinuitIF.hpp MinuitResult.hpp)

add_library(Minuit2IF ${lib_srcs} ${lib_headers})

target_link_libraries(
  Minuit2IF
  PUBLIC Core
  PRIVATE ROOT::Minuit2)

install(FILES ${lib_headers} DESTINATION include/Optimizer/Minuit2)
install(
  TARGETS Minuit2IF
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

add_subdirectory(test)
