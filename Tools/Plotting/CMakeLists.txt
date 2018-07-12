#################
# Tools Library #
#################


# add sources, headers, include directories and link libraries
# depending on the dependencies
IF( ${ROOT_FOUND} AND ${Boost_FOUND})
set(srcs ${srcs} RootPlotData.cpp)
set(headers ${headers} RootPlotData.hpp)
set(private_include_dirs ${private_include_dirs} 
	${Boost_INCLUDE_DIR}
	${ROOT_INCLUDE_DIR}
)
set(public_link_libs ${public_link_libs}
	Core
    DataReader
)
set(private_link_libs ${private_link_libs}
    ${Boost_LIBRARIES}
    ${ROOT_LIBRARIES}
)
ENDIF()


add_library (Plotting
  SHARED ${srcs}
)

target_include_directories(Plotting
  PUBLIC
    $<INSTALL_INTERFACE:/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  PRIVATE
    ${private_include_dirs}
)

target_link_libraries(Plotting
  PUBLIC
    ${public_link_libs}
  PRIVATE
    ${private_link_libs}
)

#set_target_properties(Plotting PROPERTIES LINKER_LANGUAGE CXX)

#
# Install
#
install(FILES ${headers}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Tools/Plotting    
)

install(TARGETS Plotting
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
