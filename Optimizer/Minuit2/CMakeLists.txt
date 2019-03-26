####################################
# Create Minuit2 interface library #
####################################

if( TARGET ROOT::Minuit2 )

set(lib_srcs
	MinuitFcn.cpp MinuitIF.cpp MinuitResult.cpp
)
set(lib_headers
	MinuitFcn.hpp MinuitIF.hpp MinuitResult.hpp
)

add_library(Minuit2IF
  SHARED ${lib_srcs} ${lib_headers}
)

target_link_libraries(Minuit2IF
  PUBLIC Core Boost::serialization ROOT::Minuit2
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
  message(WARNING "Target ROOT::Minuit2 is not available! Not building Minuit2 interface!")
endif()
