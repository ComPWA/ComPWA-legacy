# Create EvtGen library.

file(GLOB lib_srcs "*.cpp")
file(GLOB lib_headers "*.hh")

add_library(EvtGenIF
  SHARED ${lib_srcs}  ${lib_headers}
)

target_link_libraries(EvtGenIF
  PUBLIC
  	Core
)

target_include_directories(EvtGenIF
  PRIVATE ${QFTPP_INCLUDE_DIR}
)

install (FILES ${lib_headers} DESTINATION include/Physics/EvtGen)
install(TARGETS EvtGenIF
  LIBRARY DESTINATION lib
)
