################################
# Create Dynamics library #
################################

set(SOURCE_LIST RelativisticBreitWigner.cpp Flatte.cpp Voigtian.cpp
                Utils/Faddeeva.cc FormFactor.cpp
)

set(HEADER_LIST RelativisticBreitWigner.hpp Flatte.hpp Voigtian.hpp
                Utils/Faddeeva.hh FormFactor.hpp
)

add_library(Dynamics ${SOURCE_LIST} ${HEADER_LIST})

target_link_libraries(Dynamics PUBLIC Core FunctionTree)

target_include_directories(Dynamics PUBLIC ${Boost_INCLUDE_DIR})

#
# Install
#
install(FILES ${HEADER_LIST} DESTINATION include/Physics/Dynamics)
install(
  TARGETS Dynamics
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

#
# TESTING
#
#add_subdirectory(test)
