# This module accepts the following variables as hints for the search:
#
# ComPWA_ROOT       - The ComPWA installation prefix where to search for libraries and includes
# ComPWA_INCLUDEDIR - The directory where the ComPWA includes may be found
# ComPWA_LIBRARYDIR - The directory where the ComPWA libraries may be found
#
# This module defines the following variables as result of the search:
#
# ComPWA_INCLUDE_DIR          - The directory below which the ComPWA includes may be found
# ComPWA_LIBRARY_DIR          - The directory below which the ComPWA libraries may be found
#
# ComPWA_LIBRARIES            - A list containing the ComPWA library

# ComPWA_FOUND                - TRUE if ComPWA is found
#
#
# Example usages: FIND_PACKAGE(ComPWA) or FIND_PACKAGE(ComPWA REQUIRED)
#
# For this file the cmake version of FindBoost.cmake was used as a template and inspiration. 
###############################################################################

# Module for handling standard arguments given through FIND_PACKAGE()
INCLUDE(FindPackageHandleStandardArgs)

set(ComPWA_FOUND 0)

set(ComPWA_NAMESPACE "ComPWA") # library prefix

# If the CMake variables are not set, we check and use the corresponding
# environment variables
IF (NOT DEFINED ComPWA_ROOT AND DEFINED ENV{COMPWA_ROOT})
  SET (ComPWA_ROOT $ENV{COMPWA_ROOT})
ENDIF ()

# The more specific variables take precedence over the more generic ones
IF (ComPWA_ROOT)
  IF (NOT ComPWA_INCLUDEDIR)
    SET (ComPWA_INCLUDEDIR "${ComPWA_ROOT}/include/ComPWA")
    SET (ComPWA_INCLUDE_DIR "${ComPWA_ROOT}/include")
  ENDIF ()
  IF (NOT ComPWA_LIBRARYDIR)
    SET (ComPWA_LIBRARYDIR "${ComPWA_ROOT}/lib")
    SET (ComPWA_LIBRARY_DIR "${ComPWA_ROOT}/lib")
  ENDIF ()
ENDIF ()

# Setting FOUND based in include directory
if(ComPWA_INCLUDE_DIR)
  set(ComPWA_FOUND 1)
  # Check some files to make sure we found to correct
  # directory
endif()


set(ComPWA_INCLUDE_DIRS ${ComPWA_INCLUDE_DIR})
set(ComPWA_LIBRARY_DIRS ${ComPWA_LIBRARY_DIR})

# The above setting of ComPWA_FOUND was based only on the header files.
# Update it for the requested component libraries.
if(ComPWA_FOUND)
  # The headers were found.  Check for requested component libs.
  set(_comPWA_CHECKED_COMPONENT FALSE)
  set(_ComPWA_MISSING_COMPONENTS "")
  foreach(COMPONENT ${ComPWA_FIND_COMPONENTS})
    string(TOUPPER ${COMPONENT} UPPERCOMPONENT)

    find_library(ComPWA_${UPPERCOMPONENT}_LIBRARY NAMES ${ComPWA_NAMESPACE}_${COMPONENT}
      HINTS ${ComPWA_LIBRARY_DIRS} )

    if(ComPWA_${UPPERCOMPONENT}_LIBRARY)
      set(ComPWA_${UPPERCOMPONENT}_FOUND 1)
    else()
      set(ComPWA_${UPPERCOMPONENT}_FOUND 0)
    endif()

    set(_comPWA_CHECKED_COMPONENT TRUE)

    if(NOT ComPWA_${UPPERCOMPONENT}_FOUND AND ComPWA_FIND_REQUIRED_${COMPONENT})
      list(APPEND _ComPWA_MISSING_COMPONENTS ${COMPONENT})
    endif()
  endforeach()

  if(_ComPWA_MISSING_COMPONENTS AND _ComPWA_EXTRA_FIND_COMPONENTS)
    # Optional indirect dependencies are not counted as missing.
    list(REMOVE_ITEM _ComPWA_MISSING_COMPONENTS ${_ComPWA_EXTRA_FIND_COMPONENTS})
  endif()

  if (_ComPWA_MISSING_COMPONENTS)
    set(ComPWA_FOUND FALSE)
    # We were unable to find some libraries, so generate a sensible
    # error message that lists the libraries we were unable to find.
    string(APPEND ComPWA_ERROR_REASON "Could not find the following ComPWA libraries: ")
    foreach(COMPONENT ${_ComPWA_MISSING_COMPONENTS})
      string(TOUPPER ${COMPONENT} UPPERCOMPONENT)
      string(APPEND ComPWA_ERROR_REASON
        "        ${ComPWA_NAMESPACE}_${COMPONENT}${ComPWA_ERROR_REASON_${UPPERCOMPONENT}}")
    endforeach()

      list(LENGTH ComPWA_FIND_COMPONENTS ComPWA_NUM_COMPONENTS_WANTED)
      list(LENGTH _ComPWA_MISSING_COMPONENTS ComPWA_NUM_MISSING_COMPONENTS)
      if (${ComPWA_NUM_COMPONENTS_WANTED} EQUAL ${ComPWA_NUM_MISSING_COMPONENTS})
        string(APPEND ComPWA_ERROR_REASON
          "\n No ComPWA libraries were found. Make sure COMPWA_ROOT is set and libraries are located in $COMPWA_ROOT/lib.")
      else ()
        string(APPEND ComPWA_ERROR_REASON
          "Some (but not all) of the required ComPWA libraries were found. You may need to install these additional ComPWA libraries.")
      endif ()
  endif ()
endif()

# set(compwalibs Core DataReader RootReader Minuit2IF MinLogLH Tools qft++ HelicityFormalism DecayDynamics)
set(ComPWA_LIBRARIES)
foreach(COMPONENT ${ComPWA_FIND_COMPONENTS})
  find_library(ComPWA_${COMPONENT}_LIBRARY ComPWA_${COMPONENT} HINTS ${ComPWA_LIBRARY_DIR})
  # message( STATUS "ComPWA_${COMPONENT}_LIBRARY ComPWA_${COMPONENT} HINTS ${ComPWA_LIBRARY_DIR}")
  if(ComPWA_${COMPONENT}_LIBRARY)
    mark_as_advanced(ComPWA_${COMPONENT}_LIBRARY)
    list(APPEND ComPWA_LIBRARIES ${ComPWA_${COMPONENT}_LIBRARY})
    #list(REMOVE_ITEM ComPWA_FIND_COMPONENTS ${COMPONENT})
  endif()
endforeach()
if( ComPWA_LIBRARIES )
  list(REMOVE_DUPLICATES ComPWA_LIBRARIES)
endif()

# ------------------------------------------------------------------------
#  Add imported targets
# ------------------------------------------------------------------------

if(ComPWA_FOUND)
  foreach(COMPONENT ${ComPWA_FIND_COMPONENTS})
    if(NOT TARGET ComPWA::${COMPONENT})
      string(TOUPPER ${COMPONENT} UPPERCOMPONENT)
      if(ComPWA_${UPPERCOMPONENT}_FOUND)
        add_library(ComPWA::${COMPONENT} UNKNOWN IMPORTED)
        if(ComPWA_INCLUDE_DIRS)
          set_target_properties(ComPWA::${COMPONENT} PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${ComPWA_INCLUDE_DIRS}")
        endif()
        if(EXISTS "${ComPWA_${UPPERCOMPONENT}_LIBRARY}")
          set_target_properties(ComPWA::${COMPONENT} PROPERTIES
            IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
            IMPORTED_LOCATION "${ComPWA_${UPPERCOMPONENT}_LIBRARY}")
        endif()
      endif()
    endif()
  endforeach()
endif()

# ------------------------------------------------------------------------
#  Notification to end user about what was found
# ------------------------------------------------------------------------

set(ComPWA_LIBRARIES "")
if(ComPWA_FOUND)
  if(NOT ComPWA_FIND_QUIETLY)
    set(_comPWA_VERSION ${ComPWA_MAJOR_VERSION}.${ComPWA_MINOR_VERSION}.${ComPWA_SUBMINOR_VERSION})
    if(NOT ComPWA_FIND_COMPONENTS) # No components specified
      message(STATUS "Found ComPWA (${_comPWA_VERSION})")
    else()
      set(_component_STRING "")
      foreach( COMPONENT  ${ComPWA_FIND_COMPONENTS} )
        string( TOUPPER ${COMPONENT} UPPERCOMPONENT )
        if( ComPWA_${UPPERCOMPONENT}_FOUND )
          if(NOT ComPWA_FIND_QUIETLY)
            string(APPEND _component_STRING " ${COMPONENT}")
          endif()
          list(APPEND ComPWA_LIBRARIES ${ComPWA_${UPPERCOMPONENT}_LIBRARY})
        endif()
      endforeach()
      message(STATUS "Found ComPWA (${_comPWA_VERSION}) with the following components:\n     ${_component_STRING}" )
    endif()
  endif()

else()
  if(ComPWA_FIND_REQUIRED)
    message(SEND_ERROR "${ComPWA_ERROR_REASON}")
  else()
    if(NOT ComPWA_FIND_QUIETLY)
      # we opt not to automatically output ComPWA_ERROR_REASON here as
      # it could be quite lengthy and somewhat imposing in its requests
      # Since ComPWA is not always a required dependency we'll leave this
      # up to the end-user.
      if(ComPWA_DEBUG OR ComPWA_DETAILED_FAILURE_MSG)
        message(STATUS "Could NOT find ComPWA\n${ComPWA_ERROR_REASON}")
      else()
        message(STATUS "Could NOT find ComPWA")
      endif()
    endif()
  endif()
endif()
