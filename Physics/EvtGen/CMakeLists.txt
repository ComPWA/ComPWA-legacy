# Create EvtGen library.

file(GLOB lib_srcs "*.cpp")
file(GLOB lib_headers "*.hh")

add_library( EvtGenIF
  SHARED ${lib_srcs} ${lib_headers}
)

target_link_libraries( EvtGenIF
  Core
)

target_include_directories( EvtGenIF
  PUBLIC
)

install (FILES ${lib_headers} DESTINATION include/ComPWA)
install(TARGETS EvtGenIF
        LIBRARY DESTINATION lib
)