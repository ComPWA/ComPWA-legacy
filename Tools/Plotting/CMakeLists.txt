####################
# Plotting Library #
####################

if(${ROOT_FOUND} AND ${Boost_FOUND})

set(srcs 
  RootPlotData.cpp 
  DalitzPlot.cpp
)
set(headers
  RootPlotData.hpp
  DalitzPlot.hpp
)

add_library(Plotting
  SHARED ${srcs} ${headers}
)

target_include_directories(Plotting
  PUBLIC
    $<INSTALL_INTERFACE:/include>
  PRIVATE
    ${ROOT_INCLUDE_DIR}
)

target_link_libraries(Plotting
  PUBLIC Core Data HelicityFormalism ROOT::Core 
  PRIVATE ROOT::RIO ROOT::Tree ROOT::Physics
  	ROOT::Hist ROOT::Gpad ROOT::Graf ROOT::RIO
)

install(FILES ${headers}
    DESTINATION include/Tools/Plotting
)
install(TARGETS Plotting
    LIBRARY DESTINATION lib
)

if(${Python3_FOUND})
	install(DIRECTORY ./
  		DESTINATION lib/python/Plotting
  		FILES_MATCHING PATTERN "*.py"
  		PATTERN "__pycache__" EXCLUDE
	)
endif()

endif()
