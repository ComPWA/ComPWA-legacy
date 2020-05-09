#######################
# ComPWA Core library #
#######################

set(SOURCE_LIST
    Efficiency.cpp
    Event.cpp
    FitResult.cpp
    Logging.cpp
    ProgressBar.cpp
    Properties.cpp
    TableFormatter.cpp
    Random.cpp
)
set(HEADER_LIST
    Efficiency.hpp
    Event.hpp
    FitResult.hpp
    Logging.hpp
    FourMomentum.hpp
    ProgressBar.hpp
    Properties.hpp
    TableFormatter.hpp
    Utils.hpp
    Random.hpp
)

add_library(
  Core
  ${SOURCE_LIST}
  ${HEADER_LIST}
  $<TARGET_OBJECTS:easyloggingpp_obj>
)

target_include_directories(
  Core PUBLIC ${ELPP_INCLUDE_DIRS} Boost::serialization
)

target_link_libraries(Core PUBLIC Boost::serialization)

install(
  TARGETS Core
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Core
)

add_subdirectory(FunctionTree)
add_subdirectory(test)
