####################################
# Create Minuit2 interface library #
####################################

if( NOT ${ROOT_MINUIT2_FOUND} )
	message(WARNING "Minuit2 was not found! Not building Minuit2 interface!")
	return()
endif()

set(lib_srcs
	MinuitFcn.cpp MinuitIF.cpp MinuitResult.cpp
)
set(lib_headers
	MinuitFcn.hpp MinuitIF.hpp MinuitResult.hpp
)

add_library(Minuit2IF
  SHARED ${lib_srcs}
)

target_link_libraries(Minuit2IF
  PUBLIC
    Core
    Boost::serialization
  PRIVATE ROOT::Minuit2
)

target_include_directories( Minuit2IF
  PUBLIC ${Boost_INCLUDE_DIR}
  PRIVATE ${ROOT_INCLUDE_DIR}
)

#
# Install
#
install(FILES ${lib_headers}
    DESTINATION include/ComPWA/Optimizer/Minuit2IF
)
install(TARGETS Minuit2IF
    LIBRARY DESTINATION lib/ComPWA
)
