# Minuit2 minimizer module. Check first if library is available.
FIND_PACKAGE( Minuit2 QUIET )
IF( ${MINUIT2_FOUND} )
	ADD_SUBDIRECTORY(Minuit2)
ELSE ()
	MESSAGE(WARNING "Minuit2 not found! Not building Minuit2 interface!")
ENDIF()

# Geneva minimizer module. Check first if library is available.
FIND_PACKAGE( Geneva QUIET )
IF( ${GENEVA_FOUND} )
    ADD_SUBDIRECTORY(Geneva)
ELSE ()
	MESSAGE(WARNING "Geneva not found! Not building Geneva interface!")
ENDIF()

install (FILES Optimizer.hpp
    DESTINATION include/ComPWA
)
