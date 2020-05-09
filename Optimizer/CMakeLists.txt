set(MODULE_NAME Optimizer)

install(FILES Optimizer.hpp
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${MODULE_NAME}
)

add_subdirectory(Minuit2)

if(${GENEVA_FOUND})
  add_subdirectory(Geneva)
else()
  if(${USE_GENEVA})
    message(WARNING "Geneva not found! Not building Geneva interface!")
  else()
    message(STATUS "Flag USE_GENEVA set to OFF")
  endif()
endif()
