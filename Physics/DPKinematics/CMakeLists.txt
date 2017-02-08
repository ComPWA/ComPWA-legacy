# Create MinLogLH library.

# The REQUIRED statement is somehow not working
find_package( ROOT REQUIRED )
find_package( GSL REQUIRED gsl gslcblas )

FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)

add_library ( DPKinematics SHARED ${lib_srcs} ${lib_headers} )

#add_library ( DPKinematics SHARED ${lib_srcs} G__Dict.cxx ${lib_headers} )

target_link_libraries( DPKinematics Core ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${GSL_LIBRARIES} )
target_include_directories( DPKinematics PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} ${GSL_INCLUDE_DIR} )

### Generator dictionaries
# Currently we do not generate dictionaries of compatibility reasons
# We have to add this directory to the include directories
# so rootcint can find the header and LinkDef files
#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
#ROOT_GENERATE_DICTIONARY( G__Dict SimpleEfficiency.hpp LINKDEF libDPKinematicsLinkDef.h )
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/G__Dict.cxx.o DESTINATION lib/) 

#install (TARGETS DPKinematics DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)
