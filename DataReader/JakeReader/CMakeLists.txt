# Create JakeReader library.
find_package(ROOT COMPONENTS Core Cint RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread REQUIRED)

FILE(GLOB lib_srcs *.cpp)
add_library (JakeReader ${lib_srcs})
target_link_libraries( JakeReader ${ROOT_LIBRARIES} )

#MESSAGE( STATUS "FindROOT include dirs: " "${ROOT_INCLUDE_DIR}/" )

target_include_directories (JakeReader PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR})