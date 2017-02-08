find_package( ROOT REQUIRED )
find_package( Boost REQUIRED )

add_subdirectory( test )
add_subdirectory( Examples )

add_executable (gentwopartApp GenTwoPartApp.cpp)
target_link_libraries( gentwopartApp Core BreitWigner pthread )
target_include_directories( gentwopartApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

add_executable (dpareaTestApp DPareaTestApp.cpp)
target_link_libraries( dpareaTestApp Core DPKinematics pthread )
target_include_directories( dpareaTestApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

add_executable (genRunManager GenRunManagerApp.cpp)
target_link_libraries( genRunManager Core DPKinematics pthread )
target_include_directories( genRunManager PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

#add_executable (compareTreeAmpApp CompareTreeAmpApp.cpp)
#add_executable (lhNormalizationTestApp LHNormalizationTestApp.cpp)
#add_executable (multiVariateGaussianTestApp MultiVariateGaussianTestApp.cpp)
#add_executable (serializeApp SerializationApp.cpp)