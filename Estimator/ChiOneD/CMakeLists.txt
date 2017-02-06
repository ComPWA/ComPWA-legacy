# Create ChiOneD library.
FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)

add_library( ChiOneD SHARED ${lib_srcs} ${lib_headers} )
target_link_libraries( ChiOneD Estimator Core )

target_include_directories (ChiOneD PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

install (TARGETS ChiOneD DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)
