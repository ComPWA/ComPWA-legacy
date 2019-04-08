if(${Python3_FOUND})
  add_subdirectory(PythonInterface)

  set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
  set(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/timestamp")

  # Copy to build directory
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in ${SETUP_PY})
  add_custom_command(OUTPUT ${OUTPUT}
    COMMAND ${CMAKE_COMMAND} -E copy_directory 
      ${CMAKE_CURRENT_SOURCE_DIR}/pycompwa 
      ${CMAKE_CURRENT_BINARY_DIR}/pycompwa
    COMMAND ${CMAKE_COMMAND} -E copy
      ${CMAKE_CURRENT_SOURCE_DIR}/README.rst 
      ${CMAKE_CURRENT_BINARY_DIR}/README.rst
    COMMAND ${CMAKE_COMMAND} -E copy
      ${CMAKE_CURRENT_SOURCE_DIR}/MANIFEST.in 
      ${CMAKE_CURRENT_BINARY_DIR}/MANIFEST.in
    COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
  )
  add_custom_target(pycompwa ALL DEPENDS ${OUTPUT})

  # Copy files of the installable pycompwa module to 
  # ${CMAKE_INSTALL_PREFIX}/share/ComPWA/pycompwa. 
  install(FILES pycompwa/particle_list.xml DESTINATION share/ComPWA/)
  install(FILES MANIFEST.in README.rst DESTINATION share/ComPWA/pycompwa)
  install(FILES setup.py.in DESTINATION share/ComPWA/pycompwa RENAME setup.py)
  install(DIRECTORY pycompwa DESTINATION share/ComPWA/pycompwa)

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
