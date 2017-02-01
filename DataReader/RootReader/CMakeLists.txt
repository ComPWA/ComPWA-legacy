# Create RootReader library.
find_package(ROOT COMPONENTS Core Cint RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread REQUIRED)

FILE(GLOB lib_srcs *.cpp)
list(REMOVE_ITEM lib_srcs "${CMAKE_CURRENT_SOURCE_DIR}/DataTest.cpp")
add_library (RootReader ${lib_srcs})
target_link_libraries( RootReader ${ROOT_LIBRARIES} )

#MESSAGE( STATUS "FindROOT include dirs: " "${ROOT_INCLUDE_DIR}/" )

target_include_directories (RootReader PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR})