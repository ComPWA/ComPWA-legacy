# Find or build Minuit2
set(minuit2_mpi OFF)
set(minuit2_omp OFF)
list(
  # search paths for <package>Config.cmake files
  APPEND
  CMAKE_PREFIX_PATH
  /data/sleipnir/PANDA/PWA
  /mnt/c/install/Minuit2/latest
)
build_submodule(
  PACKAGE
  Minuit2
  CONFIG
  SUBMODULE_PATH
  ThirdParty/Minuit2
)

# Find ROOT
if(NOT ENV{ROOTSYS}) # in case thisroot.sh was not sourced
  set(ADDITIONAL_SEARCH_PATHS
      ${PROJECT_SOURCE_DIR}/ROOT
      ${PROJECT_SOURCE_DIR}/root
      /data/sleipnir/PANDA/PWA
      /mnt/c/install
  )
  find_path(
    ROOTSYS_INCLUDE TROOT.h
    HINTS ${ADDITIONAL_SEARCH_PATHS}
    PATH_SUFFIXES
      include
      latest/include
      root/include
      root/latest/include
  )
  get_filename_component(
    ROOTSYS
    ${ROOTSYS_INCLUDE}
    DIRECTORY
    CACHE
    FORCE
  )
else()
  set(ROOTSYS
      $ENV{ROOTSYS}
      CACHE FORCE
  )
endif()
list(
  # in case of ROOT config
  # https://cliutils.gitlab.io/modern-cmake/chapters/packages/ROOT.html
  APPEND
  CMAKE_PREFIX_PATH
  ${ROOTSYS}
)
find_package(
  ROOT
  COMPONENTS # for ROOT_LIBRARIES
             # see https://root.cern.ch/how/integrate-root-my-project-cmake
             Core
             EG
             Gpad
             Graf
             Hist
             MathCore
             Physics
             RIO
             Tree
  QUIET
)
if(ROOT_FOUND)
  message(
    STATUS "${BoldBlue}Found ROOT ${ROOT_VERSION} (${ROOTSYS})${ColorReset}"
  )
  get_filename_component(
    # for (target_)link_directories
    ROOT_LIBRARY_DIR
    ${ROOT_Core_LIBRARY}
    DIRECTORY
    CACHE
  )
else()
  message(
    STATUS
      "${BoldRed}ROOT not found, will not build ROOT I/O functionality${ColorReset}"
  )
endif()
