# Create ChiOneD library.
FILE(GLOB lib_srcs *.cpp)
add_library (ChiOneD SHARED ${lib_srcs})
target_link_libraries( ChiOneD Estimator ComPWACore Optimizer )

target_include_directories (ChiOneD PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

install (TARGETS ChiOneD DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)
