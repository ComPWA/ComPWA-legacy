install (FILES Optimizer.hpp
    DESTINATION include/Optimizer
)

if(TARGET ROOT::Minuit2)
  add_subdirectory(Minuit2)
else()
  message(WARNING "Target ROOT::Minuit2 is not available! Not building Minuit2 interface!")
endif()
if(${GENEVA_FOUND})
  add_subdirectory(Geneva)
else()
  message(WARNING "Geneva not found! Not building Geneva interface!")
endif()
