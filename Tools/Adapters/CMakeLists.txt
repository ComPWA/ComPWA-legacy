set(SOURCE_LIST BOSSAdapter.cpp)
set(HEADER_LIST BOSSAdapter.hpp)

add_library(Adapters ${SOURCE_LIST})

target_link_libraries(Adapters PUBLIC HelicityFormalism)

target_include_directories(Adapters PUBLIC $<INSTALL_INTERFACE:include>)

install(FILES ${HEADER_LIST}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Tools/Adapters
)
install(
  TARGETS Adapters
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
