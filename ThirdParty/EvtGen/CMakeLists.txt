# Create EvtGen library.

file(GLOB lib_srcs "*.cpp")
file(GLOB lib_headers "*.hh")

set(lib_srcs ${lib_srcs})
set(lib_headers ${lib_headers})

add_library(EvtGen
  ${lib_srcs} ${lib_headers}
)

install (FILES ${lib_headers} DESTINATION include/ThirdParty/EvtGen)
install(TARGETS EvtGen
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
