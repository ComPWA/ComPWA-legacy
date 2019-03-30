if(${Python3_FOUND})
  add_subdirectory(PythonInterface)

  set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
  set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
  set(DEPS        "${CMAKE_CURRENT_SOURCE_DIR}/pycompwa/__init__.py")
  set(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/timestamp")

  configure_file(${SETUP_PY_IN} ${SETUP_PY})

  add_custom_command(OUTPUT ${OUTPUT}
    COMMAND ${CMAKE_COMMAND} -E copy_directory 
      ${CMAKE_CURRENT_SOURCE_DIR}/pycompwa ${CMAKE_CURRENT_BINARY_DIR}/pycompwa
    COMMAND ${CMAKE_COMMAND} -E copy_directory 
      ${CMAKE_CURRENT_SOURCE_DIR}/tests ${CMAKE_CURRENT_BINARY_DIR}/tests
    COMMAND ${Python3_EXECUTABLE} ${SETUP_PY} build
    COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
    DEPENDS ${DEPS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

  add_custom_target(pycompwa ALL DEPENDS ${OUTPUT})

  # Automatic installation probably not adequate since the user could use 
  # venv, custom path etc.
  # install(CODE "execute_process(COMMAND ${Python3_EXECUTABLE} ${SETUP_PY} "
  #              "install --user WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})")

  # .. so inform the user that the python module needs to be installed separately. 
  get_filename_component(PYEXE ${Python3_EXECUTABLE} NAME)
  get_filename_component(PYSETUP ${SETUP_PY} NAME)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(Green       "${Esc}[32m")
  set(Cyan        "${Esc}[36m")
  message(STATUS "${Green}Note that the python module is build but not installed "
    "automatically. After build \n   or install you should run something like:\n "
    "${Cyan}cd ${CMAKE_CURRENT_BINARY_DIR}\n "
    "${PYEXE} ${PYSETUP} install (--user)${ColourReset}"
  )

endif()
