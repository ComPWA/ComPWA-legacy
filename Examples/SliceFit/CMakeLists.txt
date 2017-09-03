# Create executable for slice fit of J/psi -> ypipi.

find_package( ROOT REQUIRED )

FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)

add_executable ( slicefitApp SliceFitApp.cpp ${lib_headers})
target_link_libraries( slicefitApp ComPWA_AmplitudeSum ComPWA_Minuit2IF ComPWA_SliceFit ComPWA_Estimator ComPWA_RootReader ComPWA_JakeReader ComPWA_DPKinematics ComPWA_Core ${ROOT_LIBRARIES} ${Boost_LIBRARIES} pthread )
target_include_directories( slicefitApp PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

#add_executable ( slicefitrandApp SliceFitRANDApp.cpp ${lib_headers})
#target_link_libraries( slicefitrandApp ComPWA_AmplitudeSum ComPWA_Minuit2IF ComPWA_SliceFitUB ComPWA_Estimator ComPWA_RootReader ComPWA_JakeReader ComPWA_DPKinematics ComPWA_Core ${ROOT_LIBRARIES} ${Boost_LIBRARIES} pthread )
#target_include_directories( slicefitrandApp PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

add_executable ( slicefitubApp SliceFitUBApp.cpp ${lib_headers})
target_link_libraries( slicefitubApp ComPWA_AmplitudeSum ComPWA_Minuit2IF ComPWA_SliceFitUB ComPWA_Estimator ComPWA_RootReader ComPWA_JakeReader ComPWA_DPKinematics ComPWA_Core ${ROOT_LIBRARIES} ${Boost_LIBRARIES} pthread )
target_include_directories( slicefitubApp PUBLIC ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

#install (TARGETS DalitzFitD03K DESTINATION bin/ )
#install (FILES ${lib_headers} DESTINATION include/ComPWA/ )
