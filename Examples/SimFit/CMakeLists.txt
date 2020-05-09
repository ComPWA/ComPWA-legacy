if(TARGET Minuit2IF
   AND TARGET MinLogLH
   AND TARGET Tools
)

  # SimFit
  set(EXECUTABLE_NAME SimFit)
  add_executable(${EXECUTABLE_NAME} ${EXECUTABLE_NAME}.cpp)
  target_include_directories(
    ${EXECUTABLE_NAME}
    PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
            ${Boost_INCLUDE_DIR}
  )
  target_link_libraries(
    ${EXECUTABLE_NAME}
    PRIVATE Core
            MinLogLH
            Minuit2IF
            Physics
            Plotting
            RootData
  )
  install(TARGETS ${EXECUTABLE_NAME}
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
  configure_file(
    ${PROJECT_SOURCE_DIR}/Physics/particle_list.xml
    ${CMAKE_CURRENT_BINARY_DIR}/particle_list.xml
    COPYONLY
  )

endif()
