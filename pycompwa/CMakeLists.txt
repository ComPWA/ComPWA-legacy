if(${Python3_FOUND})
  add_subdirectory(PythonInterface)

  set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
  set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
  set(DEPS        "${CMAKE_CURRENT_SOURCE_DIR}/pycompwa/__init__.py")
  set(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/timestamp")

  configure_file(${SETUP_PY_IN} ${SETUP_PY})

  add_custom_command(OUTPUT ${OUTPUT}
    COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/pycompwa ${CMAKE_CURRENT_BINARY_DIR}/pycompwa
    COMMAND ${PYTHON} ${SETUP_PY} build
    COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
    DEPENDS ${DEPS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

  add_custom_target(target ALL DEPENDS ${OUTPUT})

  # install(CODE "execute_process(COMMAND ${PYTHON} ${SETUP_PY} bdist_wheel)")
  install(CODE "execute_process(COMMAND ${PYTHON} ${SETUP_PY} install --user WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})")
endif()
