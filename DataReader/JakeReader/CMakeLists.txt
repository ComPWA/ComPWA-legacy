# Create JakeReader library.
find_package(ROOT COMPONENTS Core Cint RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread REQUIRED)

FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)
add_library (ComPWA_JakeReader SHARED ${lib_srcs} ${lib_headers})
target_link_libraries( ComPWA_JakeReader ComPWA_DataReader ${ROOT_LIBRARIES} ${Boost_LIBRARIES} )

#MESSAGE( STATUS "FindROOT include dirs: " "${ROOT_INCLUDE_DIR}/" )

target_include_directories (ComPWA_JakeReader PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR})

#install (TARGETS ComPWA_JakeReader DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)
