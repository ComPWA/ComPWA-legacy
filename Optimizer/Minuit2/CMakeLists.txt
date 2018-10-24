####################################
# Create Minuit2 interface library #
####################################

if( ${ROOT_MINUIT2_FOUND} )

set(lib_srcs
	MinuitFcn.cpp MinuitIF.cpp MinuitResult.cpp
)
set(lib_headers
	MinuitFcn.hpp MinuitIF.hpp MinuitResult.hpp
)

if(${CMAKE_HOST_APPLE})
	link_directories(${ROOT_LIBRARY_DIR})
endif()

add_library(Minuit2IF
  SHARED ${lib_srcs}  ${lib_headers}
)

target_link_libraries(Minuit2IF
  PUBLIC
    Core
    Boost::serialization
    ROOT::Minuit2
)

target_include_directories( Minuit2IF
  PUBLIC ${Boost_INCLUDE_DIR} ${ROOT_INCLUDE_DIR}
)

#
# Install
#
install(FILES ${lib_headers}
  DESTINATION include/Optimizer/Minuit2
)
install(TARGETS Minuit2IF
  LIBRARY DESTINATION lib
)

else()
  message(WARNING "Minuit2 was not found! Not building Minuit2 interface!")
endif()
