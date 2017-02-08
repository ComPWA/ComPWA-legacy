# Create BreitWigner library.

# The REQUIRED statement is somehow not working
find_package( ROOT REQUIRED )

FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)
add_library ( BreitWigner SHARED ${lib_srcs} ${lib_headers})
target_link_libraries( BreitWigner Core ${ROOT_LIBRARIES} ${Boost_LIBRARIES} )
target_include_directories( BreitWigner PUBLIC ${ROOT_INCLUDE_DIR} )

#install (TARGETS BreitWigner DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)
