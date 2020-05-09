if(${ROOT_FOUND} AND ${Boost_FOUND})

  # Define library
  set(PACKAGE_NAME Plotting)
  set(SOURCE_LIST DalitzPlot.cpp RootPlotData.cpp)
  add_library(${PACKAGE_NAME} ${SOURCE_LIST})
  target_include_directories(
    ${PACKAGE_NAME}
    PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE ${ROOT_INCLUDE_DIR}
  )
  target_link_libraries(
    ${PACKAGE_NAME}
    PUBLIC Core
           Data
           Physics
           ROOT::Core
           ROOT::RIO
    PRIVATE ROOT::Gpad
            ROOT::Graf
            ROOT::Hist
            ROOT::Physics
            ROOT::Tree
  )

  # Install and export
  set(HEADER_LIST
      DalitzPlot.hpp
      HistTools.hpp
      RootPlotData.hpp
  )
  install(FILES ${HEADER_LIST}
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${MODULE_NAME}/${PACKAGE_NAME}
  )
  install(
    TARGETS ${PACKAGE_NAME}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )

endif()
