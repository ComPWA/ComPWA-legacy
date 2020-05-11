if(TARGET Minuit2IF
   AND TARGET MinLogLH
   AND TARGET Tools
   AND TARGET HelicityFormalism
)

  set(SOURCE_LIST SimFit.cpp)

  add_executable(SimFit SimFit.cpp)

  target_link_libraries(
    SimFit
    ${Boost_LIBRARIES}
    ${ROOT_LIBRARIES}
    Core
    HelicityFormalism
    MinLogLH
    Minuit2IF
    Plotting
    RootData
  )

  target_include_directories(
    SimFit
    PUBLIC $<INSTALL_INTERFACE:include>
           ${ROOT_INCLUDE_DIR}
           ${Boost_INCLUDE_DIR}
  )

  install(TARGETS SimFit RUNTIME DESTINATION bin)

  configure_file(
    ${PROJECT_SOURCE_DIR}/Physics/particle_list.xml
    ${CMAKE_CURRENT_BINARY_DIR}/particle_list.xml
    COPYONLY
  )

else()
  message(WARNING "Required targets not found! Not building\
                    SimFit executable!"
  )
endif()
