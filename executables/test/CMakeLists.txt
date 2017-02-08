# Create executable for test applications of ComPWA.

# The REQUIRED statement is somehow not working
find_package( ROOT REQUIRED )
find_package( Boost REQUIRED )

FILE(GLOB lib_srcs PolyFit.cpp)
FILE(GLOB lib_headers PolyFit.hpp)

add_library ( testApps SHARED ${lib_srcs} ${lib_headers} )
target_link_libraries( testApps ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
target_include_directories ( testApps PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

add_executable ( cx0testApp CX0TestApp.cpp ${lib_headers})
target_link_libraries( cx0testApp pthread )
target_include_directories( cx0testApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

add_executable ( treetestApp TreeTestApp.cpp ${lib_headers})
target_link_libraries( treetestApp Core pthread )
target_include_directories( treetestApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

add_executable ( dataiftestApp DataIFTestApp.cpp ${lib_headers})
target_link_libraries( dataiftestApp Core RootReader pthread )
target_include_directories( dataiftestApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

#Todo: if geneva
#add_executable ( optimizertestApp OptimizerTestApp.cpp ${lib_headers})
#target_link_libraries( optimizertestApp testApps Core RootReader MinLogLH Minuit2IF GenevaIF pthread )
#target_include_directories( optimizertestApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

add_executable ( optimizertestApp MinuitTestApp.cpp ${lib_headers})
target_link_libraries( optimizertestApp testApps Core RootReader MinLogLH Minuit2IF pthread )
target_include_directories( optimizertestApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

add_executable ( physicstestApp PhysicsTestApp.cpp ${lib_headers})
target_link_libraries( physicstestApp Core BreitWigner pthread )
target_include_directories( physicstestApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

add_executable ( estimatortestApp EstimatorTestApp.cpp ${lib_headers})
target_link_libraries( estimatortestApp Core MinLogLH ChiOneD BreitWigner RootReader pthread )
target_include_directories( estimatortestApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

add_executable ( runmanagertestApp RunManagerTestApp.cpp ${lib_headers})
target_link_libraries( runmanagertestApp Core MinLogLH BreitWigner RootReader Minuit2IF pthread )
target_include_directories( runmanagertestApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

#add_executable ( iteratorApp IteratorTestApp.cpp ${lib_headers})
#target_link_libraries( iteratorApp ${Boost_LIBRARIES} )
#target_include_directories( iteratorApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${Boost_INCLUDE_DIR})

#add_executable ( openMPTestApp OpenMPTestApp.cpp ${lib_headers})
#target_link_libraries( openMPTestApp test Core AmplitudeSum DalitzKinematics RootReader MinLogLH Minuit2IF pthread )
#target_include_directories( openMPTestApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

#add_executable ( integrationtestApp IntegrationTestApp.cpp ${lib_headers})
#target_link_libraries( integrationtestApp Core RootReader BreitWigner MinLogLH Minuit2IF pthread )
#target_include_directories( integrationtestApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

#add_executable ( integrationtestbinnedApp IntegrationTestBinnedApp.cpp ${lib_headers})
#target_link_libraries( integrationtestbinnedApp Core RootReader BreitWigner ChiOneD Minuit2IF pthread )
#target_include_directories( integrationtestbinnedApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

#install (TARGETS DalitzFitD03K DESTINATION bin/ )
#install (FILES ${lib_headers} DESTINATION include/ComPWA/ )
