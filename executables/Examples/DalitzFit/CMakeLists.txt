# Create executable for Dalitz plot fit of J/psi -> ypipi.

# The REQUIRED statement is somehow not working
find_package( ROOT REQUIRED )
find_package( Boost REQUIRED )

FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)

add_executable ( dalitzfitApp DalitzFitApp.cpp ${lib_headers})
target_link_libraries( dalitzfitApp AmplitudeSum Minuit2IF MinLogLH Estimator RootReader JakeReader DPKinematics Core ${ROOT_LIBRARIES} ${Boost_LIBRARIES} pthread )
target_include_directories( dalitzfitApp PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

add_executable ( gendalitzApp GenDalitzApp.cpp ${lib_headers})
target_link_libraries( gendalitzApp AmplitudeSum Minuit2IF MinLogLH Estimator RootReader JakeReader DPKinematics Core ${ROOT_LIBRARIES} ${Boost_LIBRARIES} pthread )
target_include_directories( gendalitzApp PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

#add_executable ( genphspApp GenPhspApp.cpp ${lib_headers})
#target_link_libraries( genphspApp AmplitudeSum Minuit2IF MinLogLH Estimator RootReader DPKinematics Core ${ROOT_LIBRARIES} ${Boost_LIBRARIES} pthread )
#target_include_directories( genphspApp PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

add_executable ( dalitzplotApp DalitzDataPlotterApp.cpp ${lib_headers})
target_link_libraries( dalitzplotApp AmplitudeSum Minuit2IF MinLogLH Estimator RootReader JakeReader DPKinematics Core ${ROOT_LIBRARIES} ${Boost_LIBRARIES} pthread )
target_include_directories( dalitzplotApp PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

#install (TARGETS dalitzfitApp gendalitzApp dalitzplotApp DESTINATION bin/ )
#install (FILES ${lib_headers} DESTINATION include/ComPWA/ )
