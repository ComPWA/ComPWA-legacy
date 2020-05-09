# Check if Git exists and if Git clone
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  set(IS_GIT_CLONE TRUE)
endif()

option(
  BUILD_ALL
  "Build all dependencies from submodules"
  OFF
)

# Add ability to print colors
include(colors)
define_colors()

# Helper function to download a Git submodule
function(add_git_submodule GIT_SUBMODULE)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive --
            ${GIT_SUBMODULE}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    RESULT_VARIABLE GIT_SUBMOD_RESULT
  )
  if(NOT
     GIT_SUBMOD_RESULT
     EQUAL
     "0"
  )
    message(
      FATAL_ERROR
        "${BoldRed}"
        "git submodule update --recursive --init ${GIT_SUBMODULE} "
        "failed with result ${GIT_SUBMOD_RESULT}, "
        "please investigate what's wrong with this Git submodule"
        "${ColorReset}"
    )
  endif()
  if(NOT
     EXISTS
     "${PROJECT_SOURCE_DIR}/${GIT_SUBMODULE}/CMakeLists.txt"
  )
    message(
      FATAL_ERROR
        "${BoldRed}"
        "Submodule ${GIT_SUBMODULE} was either not cloned correcty "
        "or does not contain a CMakeLists.txt file!\n"
        "Please update submodules and try again."
        "${ColorReset}"
    )
  endif()
  add_subdirectory(${PROJECT_SOURCE_DIR}/${GIT_SUBMODULE})
endfunction(add_git_submodule)

# Helper function for cloning and building third-party libraries
function(build_submodule)
  set(options) # function arguments
  set(oneVal PACKAGE SUBMODULE_PATH)
  set(multiVal)
  cmake_parse_arguments(
    GITSUB
    "${options}"
    "${oneVal}"
    "${multiVal}"
    ${ARGN}
  )
  option(
    BUILD_${GITSUB_PACKAGE}
    "Force building "
    ${${GITSUB_PACKAGE}}
    " from Git submodule"
  )
  find_package(${GITSUB_PACKAGE} ${GITSUB_UNPARSED_ARGUMENTS})
  if(${GITSUB_PACKAGE}_FOUND
     AND NOT BUILD_${GITSUB_PACKAGE}
     AND NOT BUILD_ALL
  )
    message(
      STATUS "${BoldBlue}"
             "Found third-party library ${GITSUB_PACKAGE}, "
             "no need to build"
             "${ColorReset}"
    )
  else()
    message(STATUS "Third-party library ${GITSUB_PACKAGE} not available "
                   "(or used option BUILD_${GITSUB_PACKAGE}=ON)"
    )
    message(
      STATUS
        "${BoldYellow}"
        "Building ${GITSUB_PACKAGE} from Git submodule ${GITSUB_SUBMODULE_PATH}"
        "${ColorReset}"
    )
    add_git_submodule(${GITSUB_SUBMODULE_PATH})
  endif()
endfunction(build_submodule)
