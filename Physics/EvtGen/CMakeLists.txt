# Create EvtGen library.
set(lib_srcs DalitzKinematics.cpp EvtGenIF.cpp)
set(lib_headers DalitzKinematics.hpp EvtGenIF.hpp)

add_library(EvtGenIF ${lib_srcs} ${lib_headers})

target_link_libraries(EvtGenIF PUBLIC Core FunctionTree HelicityFormalism
                                      EvtGen)

target_include_directories(EvtGenIF PRIVATE ${QFTPP_INCLUDE_DIR})

install(FILES ${lib_headers} DESTINATION include/Physics/EvtGen)
install(
  TARGETS EvtGenIF
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
