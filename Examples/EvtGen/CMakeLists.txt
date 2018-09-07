###########################################################
# Create executable for Dalitz plot fit of J/psi -> ypipi #
###########################################################

if(TARGET Minuit2IF AND TARGET MinLogLH
    AND TARGET Tools AND TARGET EvtGenIF
    AND ${Boost_PROGRAM_OPTIONS_FOUND}
)

set(lib_srcs EvtGen.cpp)

add_executable(EvtGen
    EvtGen.cpp)

target_link_libraries(EvtGen
    Core 
    Minuit2IF
    MinLogLH
    Tools
    EvtGenIF
    ${ROOT_LIBRARIES}
    ${Boost_LIBRARIES}
)

target_include_directories(EvtGen
  PUBLIC
    $<INSTALL_INTERFACE:include>
    ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

install(TARGETS EvtGen
  RUNTIME DESTINATION bin
)

else ()
  message( WARNING "Required targets not found! Not building\
                    EvtGen executable!")
endif()
