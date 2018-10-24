#################
# Tools Library #
#################


# add sources, headers, include directories and link libraries
# depending on the dependencies
if(${ROOT_FOUND} AND ${Boost_FOUND})

set(srcs ${srcs} ROOT/RootPlotData.cpp)
set(headers ${headers} ROOT/RootPlotData.hpp)

add_library (Plotting
  SHARED ${srcs} ${lib_headers}
)

target_include_directories(Plotting
  PUBLIC
    $<INSTALL_INTERFACE:/include>
  PRIVATE
    ${Boost_INCLUDE_DIR} ${ROOT_INCLUDE_DIR}
)

target_link_libraries(Plotting
  PUBLIC
    Core Data
  PRIVATE
    ROOT::Core ROOT::RIO ROOT::Tree ROOT::Physics
)

#
# Install
#
install(FILES ${headers}
    DESTINATION include/Tools/Plotting   
)
install(TARGETS Plotting
    LIBRARY DESTINATION lib
)

if(${Python3_FOUND})
install(DIRECTORY ./
  DESTINATION lib/python/Plotting
  FILES_MATCHING PATTERN "*.py"
  PATTERN "__pycache__" EXCLUDE
)
endif()

endif()
