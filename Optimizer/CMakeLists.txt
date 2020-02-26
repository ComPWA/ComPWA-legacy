install(FILES Optimizer.hpp DESTINATION include/Optimizer)

if(TARGET ROOT::Minuit2)
  add_subdirectory(Minuit2)
else()
  message(
    WARNING
      "Target ROOT::Minuit2 is not available! Not building Minuit2 interface!")
endif()
if(${GENEVA_FOUND})
  add_subdirectory(Geneva)
else()
  if(${USE_GENEVA})
    message(WARNING "Geneva not found! Not building Geneva interface!")
  else()
    message(STATUS "Flag USE_GENEVA set to OFF")
  endif()
endif()
