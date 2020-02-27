if(TARGET EvtGen)

  set(lib_srcs EvtGenGenerator.cpp)

  set(lib_headers EvtGenGenerator.hpp)

  add_library(EvtGenGenerator ${lib_srcs} ${lib_headers})

  target_link_libraries(EvtGenGenerator PUBLIC Core Data EvtGen
                                               HelicityFormalism)

  install(FILES ${lib_headers} DESTINATION include/Data/EvtGen)
  install(
    TARGETS EvtGenGenerator
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)

  add_subdirectory(test)

else()
  message(WARNING "EvtGen not found! Not building EvtGenGenerator!")
endif()
