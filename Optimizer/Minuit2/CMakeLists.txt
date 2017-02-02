# Create MinLogLH library.
find_package(ROOT COMPONENTS Minuit2 REQUIRED)

FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)
add_library (Minuit2 SHARED ${lib_srcs} ${lib_headers})
target_link_libraries( Minuit2 Estimator Optimizer ComPWACore ${ROOT_LIBRARIES} )
target_include_directories( Minuit2 PUBLIC ${ROOT_INCLUDE_DIR} )

install (TARGETS Minuit2 DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)
