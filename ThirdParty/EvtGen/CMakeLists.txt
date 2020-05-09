# Create EvtGen library.

file(GLOB SOURCE_LIST "*.cpp")
file(GLOB HEADER_LIST "*.hh")

set(SOURCE_LIST ${SOURCE_LIST})
set(HEADER_LIST ${HEADER_LIST})

add_library(EvtGen ${SOURCE_LIST} ${HEADER_LIST})

install(FILES ${HEADER_LIST} DESTINATION include/ThirdParty/EvtGen)
install(
  TARGETS EvtGen
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
