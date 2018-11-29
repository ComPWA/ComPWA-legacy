###########################################################
# Create executable for Dalitz plot fit of J/psi -> ypipi #
###########################################################

if(TARGET Minuit2IF AND TARGET MinLogLH
    AND TARGET Tools AND TARGET EvtGenIF AND TARGET EvtGen
    AND ${Boost_PROGRAM_OPTIONS_FOUND}
)

add_executable(EvtGenFit
    EvtGenFit.cpp)

target_link_libraries(EvtGenFit
    Core 
    Minuit2IF
    MinLogLH
    Tools
    EvtGenIF
    EvtGen
    ${ROOT_LIBRARIES}
    ${Boost_LIBRARIES}
)

target_include_directories(EvtGenFit
  PUBLIC
    $<INSTALL_INTERFACE:include>
    ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

install(TARGETS EvtGenFit
  RUNTIME DESTINATION bin
)

else ()
  message( WARNING "Required targets not found! Not building\
                    EvtGenFit executable!")
endif()
