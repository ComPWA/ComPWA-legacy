set(SOURCE_LIST MinuitIF.cpp MinuitResult.cpp)
set(HEADER_LIST
    MinuitFcn.hpp
    MinuitIF.hpp
    MinuitResult.hpp
)

add_library(Minuit2IF ${SOURCE_LIST} ${HEADER_LIST})

target_link_libraries(
  Minuit2IF
  PUBLIC Core
  PRIVATE ROOT::Minuit2
)

install(FILES ${HEADER_LIST} DESTINATION include/Optimizer/Minuit2)
install(
  TARGETS Minuit2IF
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

add_subdirectory(test)
