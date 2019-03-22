if(${Python3_FOUND})

add_test(NAME ExpertSystem-Tests
  COMMAND ${Python3_EXECUTABLE} -m pytest -v
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set_property(TEST ExpertSystem-Tests PROPERTY
  ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}")

install(DIRECTORY ./
  DESTINATION lib/python
  PATTERN "docs" EXCLUDE
  PATTERN "CMakeLists.txt" EXCLUDE
)

endif()
