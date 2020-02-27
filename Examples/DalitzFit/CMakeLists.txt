###########################################################
# Create executable for Dalitz plot fit of J/psi -> ypipi #
###########################################################

if(TARGET Minuit2IF
   AND TARGET MinLogLH
   AND TARGET Tools
   AND TARGET HelicityFormalism)

  set(lib_srcs DalitzFitApp.cpp)

  add_executable(DalitzFit DalitzFitApp.cpp)

  target_link_libraries(
    DalitzFit
    Minuit2IF
    MinLogLH
    RootData
    Integration
    Plotting
    HelicityFormalism
    Tools)

  target_include_directories(
    DalitzFit PUBLIC $<INSTALL_INTERFACE:include> ${ROOT_INCLUDE_DIR}
                     ${Boost_INCLUDE_DIR})

  install(TARGETS DalitzFit RUNTIME DESTINATION bin)

else()
  message(WARNING "Required targets not found! Not building\
                    DalitzFit executable!")
endif()
