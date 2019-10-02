if(${Python3_FOUND})

  add_subdirectory(PythonInterface)

  set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
  set(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/timestamp")

  # symlink all python package related files to the build directory
  add_custom_command(OUTPUT ${OUTPUT}
    COMMAND ${CMAKE_COMMAND} -E make_directory 
      ${CMAKE_CURRENT_BINARY_DIR}/pycompwa
    COMMAND ${CMAKE_COMMAND} -E create_symlink 
      ${CMAKE_CURRENT_SOURCE_DIR}/pycompwa/expertsystem
      ${CMAKE_CURRENT_BINARY_DIR}/pycompwa/expertsystem
    COMMAND ${CMAKE_COMMAND} -E create_symlink 
      ${CMAKE_CURRENT_SOURCE_DIR}/pycompwa/plotting
      ${CMAKE_CURRENT_BINARY_DIR}/pycompwa/plotting
    COMMAND ${CMAKE_COMMAND} -E create_symlink 
      ${CMAKE_CURRENT_SOURCE_DIR}/pycompwa/tests
      ${CMAKE_CURRENT_BINARY_DIR}/pycompwa/tests
    COMMAND ${CMAKE_COMMAND} -E create_symlink 
      ${CMAKE_CURRENT_SOURCE_DIR}/pycompwa/__init__.py
      ${CMAKE_CURRENT_BINARY_DIR}/pycompwa/__init__.py
    COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_SOURCE_DIR}/Physics/particle_list.xml
      ${CMAKE_CURRENT_BINARY_DIR}/pycompwa/particle_list.xml
    COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/README.rst 
      ${CMAKE_CURRENT_BINARY_DIR}/README.rst
    COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/MANIFEST 
      ${CMAKE_CURRENT_BINARY_DIR}/MANIFEST
    COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/setup.py
      ${SETUP_PY}
    COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
    COMMENT "Linking python files to build directory"
  )
  add_custom_target(pycompwa ALL DEPENDS ${OUTPUT})

  # Copy files of the installable pycompwa module to 
  # ${CMAKE_INSTALL_PREFIX}/share/ComPWA/pycompwa. 
  install(FILES MANIFEST README.rst DESTINATION share/ComPWA/pycompwa)
  install(FILES setup.py DESTINATION share/ComPWA/pycompwa RENAME setup.py)
  install(DIRECTORY pycompwa DESTINATION share/ComPWA/pycompwa
    PATTERN "__pycache__" EXCLUDE
    PATTERN ".pytest_cache" EXCLUDE)

  # Automatic installation of the python module probably not adequate since 
  # the user could use venv, custom path etc.
  # install(CODE "execute_process(COMMAND ${Python3_EXECUTABLE} ${SETUP_PY} "
  #              "install --user WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})")

  # ... give some hints on how to install the module
  get_filename_component(PYEXE ${Python3_EXECUTABLE} NAME)
  get_filename_component(PYSETUP ${SETUP_PY} NAME)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(Green       "${Esc}[32m")
  set(Cyan        "${Esc}[36m")
  message(STATUS "${Green}The installable python module is available in the "
    "build directory as well as in the install\ndirectory. After build or "
    "install you should run something like:\n "
    "${Cyan}cd ${CMAKE_CURRENT_BINARY_DIR} ${Green} or "
    "${Cyan}cd ${CMAKE_INSTALL_PREFIX}/share/ComPWA \n "
    "${PYEXE} ${PYSETUP} install (--user)${ColourReset}"
  )

endif()
