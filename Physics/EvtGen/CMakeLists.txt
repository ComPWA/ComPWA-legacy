# Create EvtGen library.

file(GLOB lib_srcs 
  DalitzKinematics.cpp
  EvtGenIF.cpp 
)
file(GLOB lib_headers
  DalitzKinematics.hpp
  EvtGenIF.hpp
)

set(lib_srcs ${lib_srcs} ../SubSystem.cpp ../ParticleStateTransitionKinematicsInfo.cpp)
set(lib_headers ${lib_headers} ../SubSystem.hpp)

add_library(EvtGenIF
  ${lib_srcs} ${lib_headers}
)

target_link_libraries(EvtGenIF
  PUBLIC
  	Core
  	FunctionTree
    EvtGen
)

target_include_directories(EvtGenIF
  PRIVATE 
  	${QFTPP_INCLUDE_DIR}
)

install (FILES ${lib_headers} DESTINATION include/Physics/EvtGen)
install(TARGETS EvtGenIF
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
