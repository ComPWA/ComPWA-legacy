find_package( ROOT REQUIRED )
find_package( GSL REQUIRED gslcblas gsl )

#add_executable (gentwopartApp GenTwoPartApp.cpp)
#target_link_libraries( gentwopartApp ComPWA_Core ComPWA_BreitWigner pthread ${ROOT_LIBRARIES} )
#target_include_directories( gentwopartApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR} )

#add_executable (genRunManager GenRunManagerApp.cpp)
#target_link_libraries( genRunManager ComPWA_Core ComPWA_AmplitudeSum ComPWA_RootReader ComPWA_DPKinematics pthread ${ROOT_LIBRARIES} )
#target_include_directories( genRunManager PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR})

#add_executable (compareTreeAmpApp CompareTreeAmpApp.cpp)
#target_link_libraries( compareTreeAmpApp ComPWA_Core ComPWA_RootReader ComPWA_Minuit2IF ComPWA_MinLogLH ComPWA_AmplitudeSum ${Boost_LIBRARIES} ${ROOT_LIBRARIES} )
#target_include_directories( compareTreeAmpApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${Boost_INCLUDE_DIR} ${ROOT_INCLUDE_DIR})

#add_executable (lhNormalizationTestApp LHNormalizationTestApp.cpp)
#target_link_libraries( lhNormalizationTestApp ComPWA_Core ComPWA_Minuit2IF ComPWA_RootReader ComPWA_DPKinematics ComPWA_MinLogLH ComPWA_AmplitudeSum ${Boost_LIBRARIES} ${ROOT_LIBRARIES} )
#target_include_directories( lhNormalizationTestApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${Boost_INCLUDE_DIR} ${ROOT_INCLUDE_DIR})

if( GSL_VERSION GREATER 1.9 )
#add_executable (multiVariateGaussianTestApp MultiVariateGaussianTestApp.cpp)
#target_link_libraries( multiVariateGaussianTestApp ComPWA_Core ComPWA_Minuit2IF ${ROOT_LIBRARIES} ${GSL_LIBRARIES} )
#target_include_directories( multiVariateGaussianTestApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIR} ${GSL_INCLUDE_DIR})
#MESSAGE( STATUS "GSL version " ${GSL_VERSION})

endif()

#add_executable (serializeApp SerializationApp.cpp)
#target_link_libraries( serializeApp ComPWA_Core ${Boost_LIBRARIES} )
#target_include_directories( serializeApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${Boost_INCLUDE_DIR})
