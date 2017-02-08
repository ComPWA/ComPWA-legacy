# Create Estimator library.
FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)
add_library (Estimator SHARED ${lib_srcs} ${lib_headers})

add_subdirectory (ChiOneD)
add_subdirectory (MinLogLH)
#add_subdirectory (SliceFit)
#add_subdirectory (SliceFitUB)

target_link_libraries( Estimator ${Boost_LIBRARIES} )

target_include_directories (Estimator PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

#install (TARGETS Estimator DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)
