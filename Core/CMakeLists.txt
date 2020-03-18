#######################
# ComPWA Core library #
#######################

set(lib_srcs
    Efficiency.cpp
    Event.cpp
    FitResult.cpp
    Logging.cpp
    ProgressBar.cpp
    Properties.cpp
    TableFormatter.cpp
    Random.cpp)
set(lib_headers
    Efficiency.hpp
    Event.hpp
    FitResult.hpp
    Logging.hpp
    FourMomentum.hpp
    ProgressBar.hpp
    Properties.hpp
    TableFormatter.hpp
    Utils.hpp
    Random.hpp)

add_library(Core ${lib_srcs} ${lib_headers} $<TARGET_OBJECTS:easyloggingpp_obj>)

target_include_directories(Core PUBLIC ${ELPP_INCLUDE_DIRS}
                                       Boost::serialization)

target_link_libraries(Core PUBLIC Boost::serialization)

install(
  TARGETS Core
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include/Core)

add_subdirectory(FunctionTree)
add_subdirectory(test)
