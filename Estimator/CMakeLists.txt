# Create Estimator library.
FILE(GLOB lib_srcs *.cpp)
add_library (Estimator SHARED ${lib_srcs})

add_subdirectory (ChiOneD)
add_subdirectory (MinLogLH)
#add_subdirectory (SliceFit)
#add_subdirectory (SliceFitUB)

target_link_libraries( Estimator ${Boost_LIBRARIES} ChiOneD MinLogLH )

target_include_directories (Estimator PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})