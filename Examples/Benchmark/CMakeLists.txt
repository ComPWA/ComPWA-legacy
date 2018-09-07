###########################################################
# Create executable for Dalitz plot fit of J/psi -> ypipi #
###########################################################

if( TARGET Minuit2IF AND TARGET MinLogLH
    AND TARGET Tools AND TARGET HelicityFormalism
)

set(lib_srcs Benchmark.cpp)

add_executable (Benchmark
    Benchmark.cpp)

target_link_libraries(Benchmark
    Core 
    Minuit2IF
    MinLogLH
    Tools
    HelicityFormalism
    ${ROOT_LIBRARIES}
    ${Boost_LIBRARIES}
)

target_include_directories( Benchmark
  PUBLIC 
    $<INSTALL_INTERFACE:include>
    ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

install(TARGETS Benchmark
    RUNTIME DESTINATION bin
)

else ()
  message( WARNING "Required targets not found! Not building\
                    Benchmark executable!")
endif()
