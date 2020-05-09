if(TARGET Minuit2IF
   AND TARGET MinLogLH
   AND TARGET Tools
)

  # DalitzFitApp
  set(EXECUTABLE_NAME DalitzFitApp)
  add_executable(${EXECUTABLE_NAME} ${EXECUTABLE_NAME}.cpp)
  target_include_directories(
    ${EXECUTABLE_NAME}
    PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
           ${ROOT_INCLUDE_DIR}
           ${Boost_INCLUDE_DIR}
  )
  target_link_libraries(
    ${EXECUTABLE_NAME}
    PRIVATE Integration
            MinLogLH
            Minuit2IF
            Physics
            Plotting
            RootData
            Tools
  )
  install(TARGETS ${EXECUTABLE_NAME}
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

endif()
