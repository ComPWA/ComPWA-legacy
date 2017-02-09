# Create RootReader library.
find_package( ROOT REQUIRED )
MESSAGE( STATUS "FindROOT libraries: " "${ROOT_LIBRARIES}/" )

FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)
list(REMOVE_ITEM lib_srcs "${CMAKE_CURRENT_SOURCE_DIR}/DataTest.cpp")

add_library (ComPWA_RootReader SHARED ${lib_srcs} ${lib_headers})
target_link_libraries( ComPWA_RootReader ComPWA_DataReader ${ROOT_LIBRARIES} )

target_include_directories (ComPWA_RootReader PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR})

#install (TARGETS ComPWA_RootReader DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)
