#################
# Tools Library #
#################

# Check if all requirements are found. Skip modules otherwise.
if( ${ROOT_FOUND} AND ${GSL_FOUND})

set(lib_srcs
  EvtGenGenerator.cpp
  Generate.cpp
  RootGenerator.cpp
)
set(lib_headers
  EvtGenGenerator.hpp
  FitFractions.hpp
  Generate.hpp
  ParameterTools.hpp
  PhspVolume.hpp
  RootGenerator.hpp
)

add_library(Tools
  SHARED ${lib_srcs} ${lib_headers}
)

add_library(Integration
  SHARED Integration.cpp Integration.hpp
)

target_link_libraries(Integration
  PUBLIC
  	Core
  	Data
  	pstl::ParallelSTL
)

if(${OpenMP_FOUND} AND "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  target_compile_options(Tools PRIVATE -fopenmp)
  target_compile_options(Integration PRIVATE -fopenmp)
endif()

target_link_libraries(Tools
  PUBLIC
  	Core
  	Data
  	HelicityFormalism
  	EvtGenIF
  	ROOT::Physics
  	ROOT::MathCore
  	ROOT::Core
  	pstl::ParallelSTL
    GSL::gsl
    GSL::gslcblas
  PRIVATE
  	ROOT::EG
  INTERFACE
  	${GSL_LIBRARIES}
  	ROOT::Tree
  	ROOT::Hist
  	ROOT::Gpad
  	ROOT::Graf
  	ROOT::RIO
)

target_include_directories(Tools
  PUBLIC 
    $<INSTALL_INTERFACE:include>
    ${ROOT_INCLUDE_DIR} ${GSL_INCLUDE_DIR}
)

#
# Install
#
install(FILES ${lib_headers} Integration.hpp
    DESTINATION include/Tools
)
install(TARGETS Tools Integration
    LIBRARY DESTINATION lib
)

else ()
  message(WARNING "Requirements not found! Not building Tools module!")
endif()


add_subdirectory(PythonInterface)
add_subdirectory(Plotting)
add_subdirectory(Adapters)
add_subdirectory(test)
