if(${Python3_FOUND})

install(DIRECTORY ./
  DESTINATION lib/python
  PATTERN "docs" EXCLUDE
  PATTERN "CMakeLists.txt" EXCLUDE
)

endif()