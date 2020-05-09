####################
# Plotting Library #
####################

if(${ROOT_FOUND} AND ${Boost_FOUND})

  set(SOURCE_LIST RootPlotData.cpp DalitzPlot.cpp)
  set(HEADER_LIST RootPlotData.hpp DalitzPlot.hpp)

  add_library(Plotting ${SOURCE_LIST} ${HEADER_LIST})

  target_include_directories(
    Plotting
    PUBLIC $<INSTALL_INTERFACE:/include>
    PRIVATE ${ROOT_INCLUDE_DIR}
  )

  target_link_libraries(
    Plotting
    PUBLIC Core
           Data
           HelicityFormalism
           ROOT::Core
           ROOT::RIO
    PRIVATE ROOT::Tree
            ROOT::Physics
            ROOT::Hist
            ROOT::Gpad
            ROOT::Graf
  )

  install(FILES ${HEADER_LIST} DESTINATION include/Tools/Plotting)
  install(
    TARGETS Plotting
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
  )

  if(${Python3_FOUND})
    install(
      DIRECTORY ./
      DESTINATION lib/python/Plotting
      FILES_MATCHING
      PATTERN "*.py"
      PATTERN "__pycache__" EXCLUDE
    )
  endif()

endif()
