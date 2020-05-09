set(MODULE_NAME Core)

# Define library
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
add_library(${MODULE_NAME} ${SOURCE_LIST} $<TARGET_OBJECTS:easyloggingpp>)
target_include_directories(
  ${MODULE_NAME}
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
         ${Boost_INCLUDE_DIR}
         ${ELPP_INCLUDE_DIRS}
)
target_link_libraries(${MODULE_NAME} PUBLIC Boost::serialization)

# Install and export
set(HEADER_LIST
    Efficiency.hpp
    Event.hpp
    Exceptions.hpp
    FitParameter.hpp
    FitResult.hpp
    FourMomentum.hpp
    Function.hpp
    Generator.hpp
    Kinematics.hpp
    Logging.hpp
    ProgressBar.hpp
    Properties.hpp
    Random.hpp
    Resolution.hpp
    TableFormatter.hpp
    Utils.hpp
)
install(FILES ${HEADER_LIST}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${MODULE_NAME}
)
install(
  TARGETS ${MODULE_NAME}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# Add unit tests
add_subdirectory(test)

# Add subdirectories
add_subdirectory(FunctionTree)
