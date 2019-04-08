#################
# Tools Library #
#################

# Check if all requirements are found. Skip modules otherwise.
if( ${ROOT_FOUND} AND ${GSL_FOUND})

set(lib_srcs
  EvtGenGenerator.cpp Generate.cpp RootGenerator.cpp
)
set(lib_headers
  EvtGenGenerator.hpp FitFractions.hpp Generate.hpp ParameterTools.hpp
  PhspVolume.hpp RootGenerator.hpp
)

add_library(Tools
  SHARED ${lib_srcs} ${lib_headers}
)

add_library(Integration
  SHARED Integration.cpp Integration.hpp
)

target_link_libraries(Integration
  PUBLIC Core Data pstl::ParallelSTL
)

target_link_libraries(Tools
  PUBLIC
  	Core Data	HelicityFormalism	EvtGenIF pstl::ParallelSTL 
    GSL::gsl GSL::gslcblas ROOT::MathCore ROOT::Core ROOT::Physics	
  PRIVATE
  	ROOT::EG ROOT::Tree ROOT::Hist ROOT::Gpad	ROOT::Graf ROOT::RIO
)

target_include_directories(Tools
  PUBLIC $<INSTALL_INTERFACE:include>
)

# Disable warnings. Due to the pstl dependency we see all the warnings from that.
target_compile_options(Tools PRIVATE "-w")

#
# Install
#
install(TARGETS Tools Integration
    LIBRARY DESTINATION lib
)

install(FILES ${lib_headers} Integration.hpp
    DESTINATION include/Tools
)  
else ()
  message(WARNING "Requirements not found! Not building Tools module!")
endif()


# add_subdirectory(PythonInterface)
add_subdirectory(Plotting)
add_subdirectory(Adapters)
add_subdirectory(test)
