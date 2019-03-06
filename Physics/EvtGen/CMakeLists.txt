# Create EvtGen library.

file(GLOB lib_srcs "*.cpp")
file(GLOB lib_headers "*.hpp")

set(lib_srcs ${lib_srcs} ../SubSystem.cpp ../ParticleStateTransitionKinematicsInfo.cpp)
set(lib_headers ${lib_headers} ../SubSystem.hpp)

add_library(EvtGenIF
  SHARED ${lib_srcs} ${lib_headers}
)

target_link_libraries(EvtGenIF
  PUBLIC
  	Core
  #INTERFACE
    EvtGen
)

target_include_directories(EvtGenIF
  PRIVATE 
  	${QFTPP_INCLUDE_DIR}
)

install (FILES ${lib_headers} DESTINATION include/Physics/EvtGen)
install(TARGETS EvtGenIF
  LIBRARY DESTINATION lib
)
