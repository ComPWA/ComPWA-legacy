# Create SliceFitUB library.
find_package( ROOT REQUIRED )

FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)
add_library( ComPWA_SliceFitUB SHARED ${lib_srcs} ${lib_headers} )

target_link_libraries( ComPWA_SliceFitUB ComPWA_Estimator ComPWA_Core ComPWA_DPKinematics ${ROOT_LIBRARIES} )

target_include_directories (ComPWA_SliceFitUB PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${MINUIT2}/include ${ROOT_INCLUDE_DIR})

#install (TARGETS ComPWA_SliceFitUB DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)
