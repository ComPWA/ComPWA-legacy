set(SOURCE_LIST FitFractions.cpp UpdatePTreeParameter.cpp)
set(HEADER_LIST
    FitFractions.hpp
    GoodnessOfFit.hpp
    UpdatePTreeParameter.hpp
)

add_library(Tools ${SOURCE_LIST} ${HEADER_LIST})

target_link_libraries(Tools PUBLIC Core Integration)

target_include_directories(Tools PUBLIC $<INSTALL_INTERFACE:include>)

install(
  TARGETS Tools
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

add_library(Integration Integration.cpp Integration.hpp)

target_link_libraries(
  Integration
  PUBLIC Core
         FunctionTree
         Data
         pstl::ParallelSTL
)

install(FILES ${HEADER_LIST} Integration.hpp DESTINATION include/Tools)

install(
  TARGETS Integration
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

add_subdirectory(Plotting)
add_subdirectory(Adapters)
add_subdirectory(test)
