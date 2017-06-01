###############################################################################

# Module for handling standard arguments given through FIND_PACKAGE()
INCLUDE(FindPackageHandleStandardArgs)

find_package( ROOT QUIET)

# If the CMake variables are not set, we check and use the corresponding
# environment variables

IF ((NOT DEFINED MINUIT2_ROOT) AND DEFINED ENV{MINUIT2_ROOT})
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
IF (DEFINED ROOT_INCLUDE_DIR)
	SET (PC_MINUIT2_INCLUDE_DIRS "${ROOT_INCLUDEDIR}")
ENDIF ()
IF (DEFINED ROOT_LIBRARY_DIR)
	SET (PC_MINUIT2_LIBRARY_DIRS "${ROOT_LIBRARYDIR}")
ENDIF ()

find_path(MINUIT2_INCLUDE_DIR NAMES Minuit2/FCNBase.h
   HINTS
   ${PC_MINUIT2_INCLUDEDIR}
   ${PC_MINUIT2_INCLUDE_DIRS}
   )

find_library(MINUIT2_LIBRARIES NAMES Minuit2 libMinuit2
   HINTS
   ${PC_MINUIT2_LIBRARYDIR}
   ${PC_MINUIT2_LIBRARY_DIRS}
   )


set(MINUIT2_VERSION_STRING "0.0")


FIND_PACKAGE_HANDLE_STANDARD_ARGS(Minuit2
                                  REQUIRED_VARS MINUIT2_LIBRARIES MINUIT2_INCLUDE_DIR
VERSION_VAR MINUIT2_VERSION_STRING)
