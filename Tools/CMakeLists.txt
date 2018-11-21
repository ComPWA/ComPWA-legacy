#################
# Tools Library #
#################

# Check if all requirements are found. Skip modules otherwise.
if( ${ROOT_FOUND} AND ${GSL_FOUND}
#    AND TARGET DataReader AND TARGET HelicityFormalism
)

FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)

add_library(Tools
  SHARED ${lib_srcs} ${lib_headers}
)

#target_compile_options(Tools PRIVATE -fopenmp)

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
  PRIVATE
  	ROOT::EG
  	ROOT::Tree
  	ROOT::Hist
  	ROOT::Gpad
  	ROOT::Graf
  	ROOT::RIO
  INTERFACE
  	${GSL_LIBRARIES}
)

target_include_directories(Tools
  PUBLIC 
    $<INSTALL_INTERFACE:include>
    ${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIR} ${GSL_INCLUDE_DIR}
)

#
# Install
#
install(FILES ${lib_headers}
    DESTINATION include/Tools
)
install(TARGETS Tools
    LIBRARY DESTINATION lib
)

else ()
  message(WARNING "Requirements not found! Not building Tools module!")
endif()


add_subdirectory(PythonInterface)
add_subdirectory(Plotting)
add_subdirectory(test)
