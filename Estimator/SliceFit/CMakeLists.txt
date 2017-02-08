# Create SliceFit library.
find_package( ROOT REQUIRED )

FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)
add_library( SliceFit SHARED ${lib_srcs} ${lib_headers} )

target_link_libraries( SliceFit Estimator Core ${ROOT_LIBRARIES} )

target_include_directories (SliceFit PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${MINUIT2}/include ${ROOT_INCLUDE_DIR})

#install (TARGETS SliceFit DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)
