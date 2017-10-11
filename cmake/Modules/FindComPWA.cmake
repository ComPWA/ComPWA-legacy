# This module accepts the following variables as hints for the search:
#
# COMPWA_ROOT       - The ComPWA installation prefix where to search for libraries and includes
# COMPWA_INCLUDEDIR - The directory where the ComPWA includes may be found
# COMPWA_LIBRARYDIR - The directory where the ComPWA libraries may be found
#
# This module defines the following variables as result of the search:
#
# COMPWA_INCLUDE_DIR          - The directory below which the ComPWA includes may be found
# COMPWA_LIBRARY_DIR          - The directory below which the ComPWA libraries may be found
#
# COMPWA_LIBRARIES            - A list containing the ComPWA library

# ComPWA_FOUND                - TRUE if ComPWA is found
#
#
# Example usages: FIND_PACKAGE(ComPWA) or FIND_PACKAGE(ComPWA REQUIRED)
#
###############################################################################

# Module for handling standard arguments given through FIND_PACKAGE()
INCLUDE(FindPackageHandleStandardArgs)

# If the CMake variables are not set, we check and use the corresponding
# environment variables
IF (NOT DEFINED COMPWA_ROOT AND DEFINED ENV{COMPWA_ROOT})
	SET (COMPWA_ROOT $ENV{COMPWA_ROOT})
ENDIF ()

# The more specific variables take precedence over the more generic ones
IF (COMPWA_ROOT)
	IF (NOT COMPWA_INCLUDEDIR)
	    SET (COMPWA_INCLUDEDIR "${COMPWA_ROOT}/include/ComPWA")
	    SET (COMPWA_INCLUDE_DIR "${COMPWA_ROOT}/include")
	ENDIF ()
	IF (NOT COMPWA_LIBRARYDIR)
		SET (COMPWA_LIBRARYDIR "${COMPWA_ROOT}/lib")
		SET (COMPWA_LIBRARY_DIR "${COMPWA_ROOT}/lib")
	ENDIF ()
ENDIF ()

set(compwalibs Core DataReader RootReader Minuit2IF MinLogLH Tools qft++ HelicityFormalism DecayDynamics)
set(COMPWA_LIBRARIES)
foreach(_cpt ${compwalibs} ${COMPWA_FIND_COMPONENTS})
  find_library(COMPWA_${_cpt}_LIBRARY ComPWA_${_cpt} HINTS ${COMPWA_LIBRARY_DIR})
  if(COMPWA_${_cpt}_LIBRARY)
	mark_as_advanced(COMPWA_${_cpt}_LIBRARY)
	list(APPEND COMPWA_LIBRARIES ${COMPWA_${_cpt}_LIBRARY})
	#list(REMOVE_ITEM COMPWA_FIND_COMPONENTS ${_cpt})
  endif()
endforeach()
if( COMPWA_LIBRARIES )
  list(REMOVE_DUPLICATES COMPWA_LIBRARIES)
endif()


#GET_FILENAME_COMPONENT( COMPWA_LIBRARY_DIR 
#  ${COMPWA_LIBRARIES} DIRECTORY )

# Sets ComPWA_FOUND to true if all variabels are defined
FIND_PACKAGE_HANDLE_STANDARD_ARGS( ComPWA REQUIRED_VARS
	COMPWA_LIBRARIES COMPWA_INCLUDE_DIR )

# Report some results
IF( NOT ComPWA_FIND_QUIETLY )
  IF( COMPWA_FOUND )
		MESSAGE(STATUS "ComPWA libraries: " ${COMPWA_LIBRARIES})
	ENDIF()
ENDIF()

