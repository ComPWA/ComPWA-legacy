if(${ROOT_FOUND})

  set(SOURCE_LIST RootDataIO.cpp RootEfficiency.cpp RootGenerator.cpp)

  set(HEADER_LIST RootDataIO.hpp RootEfficiency.hpp RootGenerator.hpp)

  add_library(RootData ${SOURCE_LIST} ${HEADER_LIST})

  target_link_libraries(
    RootData
    PUBLIC Core
           Data
           HelicityFormalism
           ROOT::Hist
           ROOT::MathCore
           ROOT::Core
           ROOT::Physics
    PRIVATE ROOT::EG ROOT::Tree ROOT::RIO
  )

  install(FILES ${HEADER_LIST} DESTINATION include/Data/RootData)
  install(
    TARGETS RootData
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
  )

  # Testing routines are stored in separate directory
  add_subdirectory(test)

else()
  message(WARNING "ROOT not found! Not building RootDataIO!")
endif()
