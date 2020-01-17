# - Finds ROOT instalation
# This module sets up ROOT information
# It defines:
# ROOT_FOUND             If the ROOT is found
# ROOT_INCLUDE_DIR       PATH to the include directory
# ROOT_INCLUDE_DIRS      PATH to the include directories (not cached)
# ROOT_LIBRARIES         Most common libraries
# ROOT_<name>_LIBRARY    Full path to the library <name>
# ROOT_LIBRARY_DIR       PATH to the library directory
# ROOT_ETC_DIR           PATH to the etc directory
# ROOT_DEFINITIONS       Compiler definitions
# ROOT_CXX_FLAGS         Compiler flags to used by client packages
# ROOT_C_FLAGS           Compiler flags to used by client packages
# ROOT_EXE_LINKER_FLAGS  Linker flags to used by client packages
#
# Updated by K. Smith (ksmith37@nd.edu) to properly handle
#  dependencies in ROOT_GENERATE_DICTIONARY


if("${ROOT_DIR}" STREQUAL "" AND NOT "$ENV{ROOTSYS}" STREQUAL "")
  set(ROOT_DIR "$ENV{ROOTSYS}")
endif()

if(NOT "${ROOT_DIR}" STREQUAL "")
  list(APPEND CMAKE_PREFIX_PATH ${ROOT_DIR})
endif()

find_program(ROOT_CONFIG_EXECUTABLE root-config
  HINTS ${ROOT_DIR}/bin)

execute_process(
    COMMAND ${ROOT_CONFIG_EXECUTABLE} --prefix
    OUTPUT_VARIABLE ROOTSYS
    OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(
    COMMAND ${ROOT_CONFIG_EXECUTABLE} --version
    OUTPUT_VARIABLE ROOT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(
    COMMAND ${ROOT_CONFIG_EXECUTABLE} --incdir
    OUTPUT_VARIABLE ROOT_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)
set(ROOT_INCLUDE_DIRS ${ROOT_INCLUDE_DIR})

execute_process(
    COMMAND ${ROOT_CONFIG_EXECUTABLE} --etcdir
    OUTPUT_VARIABLE ROOT_ETC_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)
set(ROOT_ETC_DIRS ${ROOT_ETC_DIR})

execute_process(
    COMMAND ${ROOT_CONFIG_EXECUTABLE} --libdir
    OUTPUT_VARIABLE ROOT_LIBRARY_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)
set(ROOT_LIBRARY_DIRS ${ROOT_LIBRARY_DIR})

set( rootlibs Core RIO EG Net Hist Graf Graf3d Gpad Tree 
     Rint Postscript Matrix Physics MathCore Thread MultiProc )

list(APPEND rootlibs ${ROOT_FIND_COMPONENTS})
set(ROOT_LIBRARIES)
foreach(_cpt ${rootlibs})
  find_library(ROOT_${_cpt}_LIBRARY ${_cpt} HINTS ${ROOT_LIBRARY_DIR})
  if(ROOT_${_cpt}_LIBRARY)
    string(TOUPPER ${_cpt} UPPERCOMPONENT)
    set(ROOT_${UPPERCOMPONENT}_FOUND TRUE)
    mark_as_advanced(ROOT_${_cpt}_LIBRARY)
    list(APPEND ROOT_LIBRARIES ${ROOT_${_cpt}_LIBRARY})
    if(ROOT_FIND_COMPONENTS)
      list(REMOVE_ITEM ROOT_FIND_COMPONENTS ${_cpt})
    endif()
  endif()
endforeach()
if(ROOT_LIBRARIES)
  list(REMOVE_DUPLICATES ROOT_LIBRARIES)
endif()

execute_process(
    COMMAND ${ROOT_CONFIG_EXECUTABLE} --cflags
    OUTPUT_VARIABLE __cflags
    OUTPUT_STRIP_TRAILING_WHITESPACE)
string(REGEX MATCHALL "-(D|U)[^ ]*" ROOT_DEFINITIONS "${__cflags}")
string(REGEX REPLACE "(^|[ ]*)-I[^ ]*" "" ROOT_CXX_FLAGS "${__cflags}")
string(REGEX REPLACE "(^|[ ]*)-I[^ ]*" "" ROOT_C_FLAGS "${__cflags}")

execute_process(
    COMMAND ${ROOT_CONFIG_EXECUTABLE} --ldflags
    OUTPUT_VARIABLE __ldflags
    OUTPUT_STRIP_TRAILING_WHITESPACE)
set(ROOT_EXE_LINKER_FLAGS "${__ldflags}")

set(ROOT_USE_FILE ${CMAKE_CURRENT_LIST_DIR}/RootUseFile.cmake)

execute_process(
  COMMAND ${ROOT_CONFIG_EXECUTABLE} --features
  OUTPUT_VARIABLE _root_options
  OUTPUT_STRIP_TRAILING_WHITESPACE)
separate_arguments(_root_options)
foreach(_opt ${_root_options})
  set(ROOT_${_opt}_FOUND TRUE)
endforeach()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ROOT DEFAULT_MSG ROOT_CONFIG_EXECUTABLE
    ROOTSYS ROOT_VERSION ROOT_INCLUDE_DIR ROOT_LIBRARIES ROOT_LIBRARY_DIR)

mark_as_advanced(ROOT_CONFIG_EXECUTABLE)

include(CMakeParseArguments)
find_program(ROOTCLING_EXECUTABLE rootcling HINTS ${ROOT_DIR}/bin)
find_program(GENREFLEX_EXECUTABLE genreflex HINTS ${ROOT_DIR}/bin)
find_package(GCCXML)


# ------------------------------------------------------------------------
#  Add imported targets
# ------------------------------------------------------------------------

set(_ROOT_IMPORTED_TARGETS TRUE)
if(ROOT_FOUND)
  foreach(COMPONENT ${rootlibs} ${ROOT_FIND_COMPONENTS})
    if(_ROOT_IMPORTED_TARGETS AND NOT TARGET ROOT::${COMPONENT})
      string(TOUPPER ${COMPONENT} UPPERCOMPONENT)
      if(ROOT_${UPPERCOMPONENT}_FOUND)
				# message(STATUS "  ${COMPONENT}")
        set(_FOUND_COMPONENTS "${_FOUND_COMPONENTS} ${COMPONENT}")
        if(ROOT_USE_STATIC_LIBS)
          add_library(ROOT::${COMPONENT} STATIC IMPORTED)
        else()
          # Even if ROOT_USE_STATIC_LIBS is OFF, we might have static
          # libraries as a result.
          add_library(ROOT::${COMPONENT} UNKNOWN IMPORTED)
        endif()
        if(ROOT_INCLUDE_DIRS)
          set_target_properties(ROOT::${COMPONENT} PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${ROOT_INCLUDE_DIRS}")
        endif()
        if(EXISTS "${ROOT_${COMPONENT}_LIBRARY}")
          set_target_properties(ROOT::${COMPONENT} PROPERTIES
            IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
            IMPORTED_LOCATION "${ROOT_${COMPONENT}_LIBRARY}")
        endif()
        if(EXISTS "${ROOT_${COMPONENT}_LIBRARY_RELEASE}")
          set_property(TARGET ROOT::${COMPONENT} APPEND PROPERTY
            IMPORTED_CONFIGURATIONS RELEASE)
          set_target_properties(ROOT::${COMPONENT} PROPERTIES
            IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
            IMPORTED_LOCATION_RELEASE "${ROOT_${COMPONENT}_LIBRARY}")
        endif()
	      if(EXISTS "${ROOT_${COMPONENT}_LIBRARY_DEBUG}")
          set_property(TARGET ROOT::${COMPONENT} APPEND PROPERTY
            IMPORTED_CONFIGURATIONS DEBUG)
          set_target_properties(ROOT::${COMPONENT} PROPERTIES
            IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
            IMPORTED_LOCATION_DEBUG "${ROOT_${COMPONENT}_LIBRARY}")
        endif()
        if(_ROOT_${UPPERCOMPONENT}_DEPENDENCIES)
          unset(_ROOT_${UPPERCOMPONENT}_TARGET_DEPENDENCIES)
          foreach(dep ${_ROOT_${UPPERCOMPONENT}_DEPENDENCIES})
            list(APPEND _ROOT_${UPPERCOMPONENT}_TARGET_DEPENDENCIES ROOT::${dep})
          endforeach()
          if(COMPONENT STREQUAL "thread")
            list(APPEND _ROOT_${UPPERCOMPONENT}_TARGET_DEPENDENCIES Threads::Threads)
          endif()
          set_target_properties(ROOT::${COMPONENT} PROPERTIES
            INTERFACE_LINK_LIBRARIES "${_ROOT_${UPPERCOMPONENT}_TARGET_DEPENDENCIES}")
        endif()
        if(_ROOT_${UPPERCOMPONENT}_COMPILER_FEATURES)
          set_target_properties(ROOT::${COMPONENT} PROPERTIES
            INTERFACE_COMPILE_FEATURES "${_ROOT_${UPPERCOMPONENT}_COMPILER_FEATURES}")
        endif()
      endif()
    endif()
  endforeach()
  message(STATUS "ROOT version: ${ROOT_VERSION} at ${ROOT_LIBRARY_DIRS}" )
  message(STATUS "Found the following ROOT libraries:" )
  message(STATUS " ${_FOUND_COMPONENTS}" )
endif()


#----------------------------------------------------------------------------
# function ROOT_GENERATE_DICTIONARY( dictionary
#                                    header1 header2 ...
#                                    LINKDEF linkdef1 ...
#                                    OPTIONS opt1...)
function(ROOT_GENERATE_DICTIONARY dictionary)
  CMAKE_PARSE_ARGUMENTS(ARG "" "" "LINKDEF;OPTIONS" "" ${ARGN})
  #---Get the list of include directories------------------
  get_directory_property(incdirs INCLUDE_DIRECTORIES)
  set(includedirs)
  foreach( d ${incdirs})
     set(includedirs ${includedirs} -I${d})
  endforeach()
  #---Get the list of header files-------------------------
  set(headerfiles)
  foreach(fp ${ARG_UNPARSED_ARGUMENTS})
    if(${fp} MATCHES "[*?]") # Is this header a globbing expression?
      file(GLOB files ${fp})
      foreach(f ${files})
        if(NOT f MATCHES LinkDef) # skip LinkDefs from globbing result
          set(headerfiles ${headerfiles} ${f})
        endif()
      endforeach()
    else()
      find_file(headerFile ${fp} HINTS ${incdirs})
      set(headerfiles ${headerfiles} ${headerFile})
      unset(headerFile CACHE)
    endif()
  endforeach()
  #---Get LinkDef.h file------------------------------------
  set(linkdefs)
  foreach( f ${ARG_LINKDEF})
    find_file(linkFile ${f} HINTS ${incdirs} "./")
    set(linkdefs ${linkdefs} ${linkFile})
    unset(linkFile CACHE)
  endforeach()
  #---call rootcling------------------------------------------
  add_custom_command(OUTPUT ${dictionary}.cxx
                     COMMAND rootcling -f ${dictionary}.cxx 
		     ${ARG_OPTIONS} ${includedirs} ${headerfiles} ${linkdefs}
                     DEPENDS ${headerfiles} ${linkdefs} 
		     COMMENT "Generating root dictionary ${dictionary}"
		     )
endfunction()

#----------------------------------------------------------------------------
# function REFLEX_GENERATE_DICTIONARY(dictionary
#                                     header1 header2 ...
#                                     SELECTION selectionfile ...
#                                     OPTIONS opt1...)
function(REFLEX_GENERATE_DICTIONARY dictionary)
  CMAKE_PARSE_ARGUMENTS(ARG "" "" "SELECTION;OPTIONS" "" ${ARGN})
  #---Get the list of header files-------------------------
  set(headerfiles)
  foreach(fp ${ARG_UNPARSED_ARGUMENTS})
    file(GLOB files ${fp})
    if(files)
      foreach(f ${files})
        set(headerfiles ${headerfiles} ${f})
      endforeach()
    else()
      set(headerfiles ${headerfiles} ${fp})
    endif()
  endforeach()
  #---Get Selection file------------------------------------
  if(IS_ABSOLUTE ${ARG_SELECTION})
    set(selectionfile ${ARG_SELECTION})
  else()
    set(selectionfile ${CMAKE_CURRENT_SOURCE_DIR}/${ARG_SELECTION})
  endif()
  #---Get the list of include directories------------------
  get_directory_property(incdirs INCLUDE_DIRECTORIES)
  set(includedirs)
  foreach( d ${incdirs})
    set(includedirs ${includedirs} -I${d})
  endforeach()
  #---Get preprocessor definitions--------------------------
  get_directory_property(defs COMPILE_DEFINITIONS)
  foreach( d ${defs})
   set(definitions ${definitions} -D${d})
  endforeach()
  #---Nanes and others---------------------------------------
  set(gensrcdict ${dictionary}.cpp)
  if(MSVC)
    set(gccxmlopts "--gccxmlopt=\"--gccxml-compiler cl\"")
  else()
    #set(gccxmlopts "--gccxmlopt=\'--gccxml-cxxflags -m64 \'")
    set(gccxmlopts)
  endif()
  #set(rootmapname ${dictionary}Dict.rootmap)
  #set(rootmapopts --rootmap=${rootmapname} --rootmap-lib=${libprefix}${dictionary}Dict)
  #---Check GCCXML and get path-----------------------------
  if(GCCXML)
    get_filename_component(gccxmlpath ${GCCXML} PATH)
  else()
    message(WARNING "GCCXML not found. Install and setup your environment to find 'gccxml' executable")
  endif()
  #---Actual command----------------------------------------
  add_custom_command(OUTPUT ${gensrcdict} ${rootmapname}
                     COMMAND ${GENREFLEX_EXECUTABLE} ${headerfiles} -o ${gensrcdict} ${gccxmlopts} ${rootmapopts} --select=${selectionfile}
                             --gccxmlpath=${gccxmlpath} ${ARG_OPTIONS} ${includedirs} ${definitions}
                     DEPENDS ${headerfiles} ${selectionfile})
endfunction()

