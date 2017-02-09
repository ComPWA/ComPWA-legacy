# Create BreitWigner library.

# The REQUIRED statement is somehow not working
find_package( ROOT REQUIRED )

FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)
add_library ( ComPWA_BreitWigner SHARED ${lib_srcs} ${lib_headers})
target_link_libraries( ComPWA_BreitWigner ComPWA_Core ${ROOT_LIBRARIES} ${Boost_LIBRARIES} )
target_include_directories( ComPWA_BreitWigner PUBLIC ${ROOT_INCLUDE_DIR} )

#install (TARGETS BreitWigner DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)
