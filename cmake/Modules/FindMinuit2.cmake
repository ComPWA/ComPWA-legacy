# This module accepts the following variables as hints for the search:
#
# MINUIT2_ROOT       - The Minuit2 installation prefix where to search for libraries and includes
# MINUIT2_INCLUDEDIR - The directory where the Minuit2 includes may be found
# MINUIT2_LIBRARYDIR - The directory where the Minuit2 libraries may be found
#
# This module defines the following variables as result of the search:
#
# MINUIT2_INCLUDE_DIR          - The directory below which the Minuit2 includes may be found
# MINUIT2_LIBRARY_DIR          - The directory below which the Minuit2 libraries may be found
#
# MINUIT2_LIBRARIES            - A list containing the Minuit2 library

# Minuit2_FOUND                - TRUE if Minuit2 is found
#
#
# Example usages: FIND_PACKAGE(Minuit2) or FIND_PACKAGE(Minuit2 REQUIRED)
#
###############################################################################

# Module for handling standard arguments given through FIND_PACKAGE()
INCLUDE(FindPackageHandleStandardArgs)

# If the CMake variables are not set, we check and use the corresponding
# environment variables
IF (NOT DEFINED MINUIT2_ROOT AND DEFINED ENV{MINUIT2_ROOT})
	SET (MINUIT2_ROOT $ENV{MINUIT2_ROOT})
ENDIF ()

# The more specific variables take precedence over the more generic ones
IF (MINUIT2_ROOT)
	IF (NOT MINUIT2_INCLUDEDIR)
		SET (PC_MINUIT2_INCLUDEDIR "${MINUIT2_ROOT}/include")
	ENDIF ()
	IF (NOT MINUIT2_LIBRARYDIR)
		SET (PC_MINUIT2_LIBRARYDIR "${MINUIT2_ROOT}/lib")
	ENDIF ()
ENDIF ()

# Try to find Minuit2 within the ROOT installation
find_package( ROOT QUIET )
IF( ROOT_FOUND )
	IF (DEFINED ROOT_INCLUDE_DIR)
		IF (NOT PC_MINUIT2_INCLUDEDIR )
			SET (PC_MINUIT2_INCLUDE_DIRS "${ROOT_INCLUDE_DIR}")
		ENDIF ()
	ENDIF ()
	IF (DEFINED ROOT_LIBRARY_DIR)
		SET (PC_MINUIT2_LIBRARY_DIRS "${ROOT_LIBRARY_DIR}")
	ENDIF ()
ENDIF()

FIND_PATH(MINUIT2_INCLUDE_DIR NAMES Minuit2/FCNBase.h
	HINTS
	${PC_MINUIT2_INCLUDEDIR}
	${PC_MINUIT2_INCLUDE_DIRS}
	)

FIND_LIBRARY(MINUIT2_LIBRARIES NAMES Minuit2
	HINTS
	${PC_MINUIT2_LIBRARYDIR}
	${PC_MINUIT2_LIBRARY_DIRS}
	)

GET_FILENAME_COMPONENT( MINUIT2_LIBRARY_DIR 
  ${MINUIT2_LIBRARIES} DIRECTORY )

# Sets Minuit2_FOUND to true if all variabels are defined
FIND_PACKAGE_HANDLE_STANDARD_ARGS( Minuit2 REQUIRED_VARS
	MINUIT2_LIBRARIES MINUIT2_INCLUDE_DIR )

# Report some results
IF( NOT Minuit2_FIND_QUIETLY )
  IF( MINUIT2_FOUND )
		MESSAGE(STATUS "Minuit2 libraries: " ${MINUIT2_LIBRARIES})
	#ELSE()
		#MESSAGE(STATUS "Minuit2 installation could not be found!")
	ENDIF()
ENDIF()

