# Minuit2 minimizer module. Check first if library is available.
FIND_PACKAGE( Minuit2 QUIET )
MESSAGE(STATUS "ddd " ${Minuit2_FOUND} "--" ${MINUIT2_LIBRARIES})
IF( ${Minuit2_FOUND} )
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

MESSAGE(STATUS " -- " ${Minuit2_FOUND} " -- " ${GENEVA_FOUND})
