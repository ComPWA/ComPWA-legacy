set(lib_srcs BOSSAdapter.cpp)
set(lib_headers BOSSAdapter.hpp)

add_library(Adapters ${lib_srcs})

target_link_libraries(Adapters PUBLIC HelicityFormalism)

target_include_directories(Adapters PUBLIC $<INSTALL_INTERFACE:include>)

install(FILES ${lib_headers} DESTINATION include/Tools/Adapters)
install(
  TARGETS Adapters
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
