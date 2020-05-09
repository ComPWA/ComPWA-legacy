if(ROOT_FOUND)

  # Define library
  set(PACKAGE_NAME RootData)
  set(SOURCE_LIST
      RootDataIO.cpp
      RootEfficiency.cpp
      RootGenerator.cpp
  )
  add_library(${PACKAGE_NAME} ${SOURCE_LIST})
  target_include_directories(
    ${PACKAGE_NAME}
    PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
           ${ROOT_INCLUDE_DIR}
  )
  target_link_libraries(
    ${PACKAGE_NAME}
    PUBLIC Core
           Data
           Physics
           ROOT::Hist
           ROOT::MathCore
           ROOT::Core
           ROOT::Physics
    PRIVATE ROOT::EG
            ROOT::Tree
            ROOT::RIO
  )

  # Install and export
  set(HEADER_LIST
      RootDataIO.hpp
      RootEfficiency.hpp
      RootGenerator.hpp
  )
  install(FILES ${HEADER_LIST}
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${MODULE_NAME}/Root
  )
  install(
    TARGETS ${PACKAGE_NAME}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )

  # Add unit tests
  add_subdirectory(test)

endif()
