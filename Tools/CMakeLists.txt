#################
# Tools Library #
#################

# Check if all requirements are found. Skip modules otherwise.
if(${GSL_FOUND})

set(lib_srcs
  UpdatePTreeParameter.cpp
)
set(lib_headers
  FitFractions.hpp
  PhspVolume.hpp 
  UpdatePTreeParameter.hpp
)

add_library(Tools
  ${lib_srcs} ${lib_headers}
)

target_link_libraries(Tools
  PUBLIC
  	Core Integration
    GSL::gsl GSL::gslcblas
)

target_include_directories(Tools
  PUBLIC $<INSTALL_INTERFACE:include>
)

install(TARGETS Tools
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

else ()
  message(WARNING "Requirements not found! Not building Tools module!")
endif()

add_library(Integration
  Integration.cpp Integration.hpp
)

target_link_libraries(Integration
  PUBLIC Core FunctionTree Data pstl::ParallelSTL
)

install(FILES ${lib_headers} Integration.hpp
    DESTINATION include/Tools
)  

install(TARGETS Integration
    LIBRARY DESTINATION lib
)

add_subdirectory(Plotting)
add_subdirectory(Adapters)
add_subdirectory(test)
