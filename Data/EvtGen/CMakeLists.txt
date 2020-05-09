if(TARGET EvtGen)

  set(SOURCE_LIST EvtGenGenerator.cpp)

  set(HEADER_LIST EvtGenGenerator.hpp)

  add_library(EvtGenGenerator ${SOURCE_LIST} ${HEADER_LIST})

  target_link_libraries(
    EvtGenGenerator
    PUBLIC Core
           Data
           EvtGen
           HelicityFormalism
  )

  install(FILES ${HEADER_LIST}
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Data/EvtGen
  )
  install(
    TARGETS EvtGenGenerator
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )

  add_subdirectory(test)

else()
  message(WARNING "EvtGen not found! Not building EvtGenGenerator!")
endif()
