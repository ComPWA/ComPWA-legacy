# Create ChiOneD library.
FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)

add_library( ComPWA_ChiOneD SHARED ${lib_srcs} ../Estimator.hpp ${lib_headers} )
target_link_libraries( ComPWA_ChiOneD ComPWA_Estimator ComPWA_Core )

target_include_directories (ComPWA_ChiOneD PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

#install (TARGETS ComPWA_ChiOneD DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)
