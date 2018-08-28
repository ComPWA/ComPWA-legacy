#################
# Tools Library #
#################


# add sources, headers, include directories and link libraries
# depending on the dependencies
if(${ROOT_FOUND} AND ${Boost_FOUND})

set(srcs ${srcs} ROOT/RootPlotData.cpp)
set(headers ${headers} ROOT/RootPlotData.hpp)

add_library (Plotting
  SHARED ${srcs}
)

target_include_directories(Plotting
  PUBLIC
    $<INSTALL_INTERFACE:/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  PRIVATE
    ${Boost_INCLUDE_DIR} ${ROOT_INCLUDE_DIR}
)

target_link_libraries(Plotting
  PUBLIC
    Core DataReader
  PRIVATE
    ROOT::Core ROOT::RIO ROOT::Tree ROOT::Physics
)

#
# Install
#
install(FILES ${headers}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Tools/Plotting   
)
install(TARGETS Plotting
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}    
)

ENDIF()
