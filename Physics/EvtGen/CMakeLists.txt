# Create EvtGen library.
set(SOURCE_LIST DalitzKinematics.cpp EvtGenIF.cpp)
set(HEADER_LIST DalitzKinematics.hpp EvtGenIF.hpp)

add_library(EvtGenIF ${SOURCE_LIST} ${HEADER_LIST})

target_link_libraries(
  EvtGenIF
  PUBLIC Core
         FunctionTree
         HelicityFormalism
         EvtGen
)

target_include_directories(EvtGenIF PRIVATE ${QFTPP_INCLUDE_DIR})

install(FILES ${HEADER_LIST}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Physics/EvtGen
)
install(
  TARGETS EvtGenIF
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
