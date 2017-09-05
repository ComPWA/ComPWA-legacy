# Minuit2 minimizer module. Check first if 
# library is available.
find_package( Minuit2 )
IF( DEFINED Minuit2_FOUND )
	add_subdirectory(Minuit2)
ELSE ()
	message(WARNING "Minuit2 not found! Not building Minuit2 interface!")
ENDIF()

# Geneva minimizer module. Check first if 
# library is available.
find_package( Geneva QUIET )
IF( DEFINED GENEVA_FOUND )
    add_subdirectory(Geneva)
ELSE ()
	message(WARNING "Geneva not found! Not building Geneva interface!")
ENDIF()
