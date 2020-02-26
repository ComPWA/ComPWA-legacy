if(${ROOT_FOUND})

  set(lib_srcs RootDataIO.cpp RootEfficiency.cpp RootGenerator.cpp)

  set(lib_headers RootDataIO.hpp RootEfficiency.hpp RootGenerator.hpp)

  add_library(RootData ${lib_srcs} ${lib_headers})

  target_link_libraries(
    RootData
    PUBLIC Core
           Data
           HelicityFormalism
           ROOT::Hist
           ROOT::MathCore
           ROOT::Core
           ROOT::Physics
    PRIVATE ROOT::EG ROOT::Tree ROOT::RIO)

  install(FILES ${lib_headers} DESTINATION include/Data/RootData)
  install(
    TARGETS RootData
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)

  # Testing routines are stored in separate directory
  add_subdirectory(test)

else()
  message(WARNING "ROOT not found! Not building RootDataIO!")
endif()
