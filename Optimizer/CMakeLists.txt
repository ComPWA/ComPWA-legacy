# Create Optimizer libraries.
#FILE(GLOB lib_srcs *.cpp)
#FILE(GLOB lib_headers *.hpp)
#add_library (ComPWA_Optimizer SHARED ${lib_srcs} ${lib_headers})
find_package( Geneva QUIET )

add_subdirectory(Minuit2)
IF(GENEVA_FOUND)
    add_subdirectory(Geneva)
ENDIF()

#target_link_libraries( ComPWA_Optimizer ComPWA_Core ${Boost_LIBRARIES} )

#target_include_directories(ComPWA_Optimizer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

#install (TARGETS ComPWA_Optimizer DESTINATION lib/ComPWA)
#install (FILES ${lib_headers} DESTINATION include/ComPWA)
