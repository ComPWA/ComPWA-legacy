# Create executable for test applications of ComPWA.

# The REQUIRED statement is somehow not working
find_package( ROOT REQUIRED )
find_package( Boost REQUIRED )
#find_package( Geneva )

FILE(GLOB lib_srcs PolyFit.cpp)
FILE(GLOB lib_headers PolyFit.hpp)

add_library ( ComPWA_testApps SHARED ${lib_srcs} ${lib_headers} )
target_link_libraries( ComPWA_testApps ComPWA_Core ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
target_include_directories ( ComPWA_testApps PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

add_executable ( cx0testApp CX0TestApp.cpp ${lib_headers})
target_link_libraries( cx0testApp pthread )
target_include_directories( cx0testApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

add_executable ( treetestApp TreeTestApp.cpp ${lib_headers})
target_link_libraries( treetestApp ComPWA_Core pthread )
target_include_directories( treetestApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

add_executable ( dataiftestApp DataIFTestApp.cpp ${lib_headers})
target_link_libraries( dataiftestApp ComPWA_Core ComPWA_RootReader ${ROOT_LIBRARIES} pthread )
target_include_directories( dataiftestApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR} )

IF(GENEVA_FOUND)
    add_executable ( optimizertestApp OptimizerTestApp.cpp ${lib_headers})
    target_link_libraries( optimizertestApp ComPWA_testApps ComPWA_Core ComPWA_RootReader ComPWA_MinLogLH ComPWA_Minuit2IF ComPWA_GenevaIF ${Boost_LIBRARIES} pthread )
    target_include_directories( optimizertestApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
ENDIF()

add_executable ( minuittestApp MinuitTestApp.cpp ${lib_headers})
target_link_libraries( minuittestApp ComPWA_testApps ComPWA_Core ComPWA_RootReader ComPWA_MinLogLH ComPWA_Minuit2IF ${Boost_LIBRARIES} pthread )
target_include_directories( minuittestApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

add_executable ( physicstestApp PhysicsTestApp.cpp ${lib_headers})
target_link_libraries( physicstestApp ComPWA_Core ComPWA_BreitWigner ${ROOT_LIBRARIES} pthread )
target_include_directories( physicstestApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR} )

add_executable ( estimatortestApp EstimatorTestApp.cpp ${lib_headers})
target_link_libraries( estimatortestApp ComPWA_Core ComPWA_MinLogLH ComPWA_ChiOneD ComPWA_BreitWigner ComPWA_RootReader ${Boost_LIBRARIES} pthread )
target_include_directories( estimatortestApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

add_executable ( runmanagertestApp RunManagerTestApp.cpp ${lib_headers})
target_link_libraries( runmanagertestApp ComPWA_Core ComPWA_MinLogLH ComPWA_BreitWigner ComPWA_RootReader ComPWA_Minuit2IF ${Boost_LIBRARIES} pthread )
target_include_directories( runmanagertestApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

#add_executable ( iteratorApp IteratorTestApp.cpp ${lib_headers})
#target_link_libraries( iteratorApp ${Boost_LIBRARIES} )
#target_include_directories( iteratorApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${Boost_INCLUDE_DIR})

#add_executable ( openMPTestApp OpenMPTestApp.cpp ${lib_headers})
#target_link_libraries( openMPTestApp test ComPWA_Core ComPWA_AmplitudeSum ComPWA_DalitzKinematics ComPWA_RootReader ComPWA_MinLogLH ComPWA_Minuit2IF ${Boost_LIBRARIES} pthread )
#target_include_directories( openMPTestApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

#add_executable ( integrationtestApp IntegrationTestApp.cpp ${lib_headers})
#target_link_libraries( integrationtestApp ComPWA_Core ComPWA_RootReader ComPWA_BreitWigner ComPWA_MinLogLH ComPWA_Minuit2IF ${Boost_LIBRARIES} pthread )
#target_include_directories( integrationtestApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

#add_executable ( integrationtestbinnedApp IntegrationTestBinnedApp.cpp ${lib_headers})
#target_link_libraries( integrationtestbinnedApp ComPWA_Core ComPWA_RootReader ComPWA_BreitWigner ComPWA_ChiOneD MComPWA_inuit2IF ${Boost_LIBRARIES} pthread )
#target_include_directories( integrationtestbinnedApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

#install (TARGETS DalitzFitD03K DESTINATION bin/ )
#install (FILES ${lib_headers} DESTINATION include/ComPWA/ )
