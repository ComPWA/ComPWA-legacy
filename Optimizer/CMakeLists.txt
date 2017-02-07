# Create Optimizer library.
#FILE(GLOB lib_srcs *.cpp)
#FILE(GLOB lib_headers *.hpp)
#add_library (Optimizer SHARED ${lib_srcs} ${lib_headers})

add_subdirectory(Minuit2)
#add_subdirectory(Geneva)

#target_link_libraries( Optimizer Core ${Boost_LIBRARIES} )

#target_include_directories(Optimizer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

#install (TARGETS Optimizer DESTINATION lib/ComPWA)
#install (FILES ${lib_headers} DESTINATION include/ComPWA)
