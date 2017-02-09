# Create executable for Dalitz plot fit of J/psi -> ypipi.

# The REQUIRED statement is somehow not working
find_package( ROOT REQUIRED )
find_package( Boost REQUIRED )

FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)

add_executable ( dalitzfitApp DalitzFitApp.cpp ${lib_headers})
target_link_libraries( dalitzfitApp ComPWA_AmplitudeSum ComPWA_Minuit2IF ComPWA_MinLogLH ComPWA_Estimator ComPWA_RootReader ComPWA_JakeReader ComPWA_DPKinematics ComPWA_Core ${ROOT_LIBRARIES} ${Boost_LIBRARIES} pthread )
target_include_directories( dalitzfitApp PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

add_executable ( gendalitzApp GenDalitzApp.cpp ${lib_headers})
target_link_libraries( gendalitzApp ComPWA_AmplitudeSum ComPWA_Minuit2IF ComPWA_MinLogLH ComPWA_Estimator ComPWA_RootReader ComPWA_JakeReader ComPWA_DPKinematics ComPWA_Core ${ROOT_LIBRARIES} ${Boost_LIBRARIES} pthread )
target_include_directories( gendalitzApp PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

#add_executable ( genphspApp GenPhspApp.cpp ${lib_headers})
#target_link_libraries( genphspApp ComPWA_AmplitudeSum ComPWA_Minuit2IF ComPWA_MinLogLH ComPWA_Estimator ComPWA_RootReader ComPWA_DPKinematics ComPWA_Core ${ROOT_LIBRARIES} ${Boost_LIBRARIES} pthread )
#target_include_directories( genphspApp PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

add_executable ( dalitzplotApp DalitzDataPlotterApp.cpp ${lib_headers})
target_link_libraries( dalitzplotApp ComPWA_AmplitudeSum ComPWA_Minuit2IF ComPWA_MinLogLH ComPWA_Estimator ComPWA_RootReader ComPWA_JakeReader ComPWA_DPKinematics ComPWA_Core ${ROOT_LIBRARIES} ${Boost_LIBRARIES} pthread )
target_include_directories( dalitzplotApp PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

#install (TARGETS dalitzfitApp gendalitzApp dalitzplotApp DESTINATION bin/ )
#install (FILES ${lib_headers} DESTINATION include/ComPWA/ )
