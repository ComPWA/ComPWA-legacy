#
# Easyloggingpp - Logging library
#

# Create object file
add_library(easyloggingpp SHARED easylogging++.cc)

#target_compile_options( easyloggingpp_object PRIVATE -w )

#target_include_directories(easyloggingpp_object PRIVATE
#                           ${CMAKE_CURRENT_SOURCE_DIR})

# Add ${EASYLOGGINGPP_ENABLE} to TARGET_LINK_LIBRARIES for earch
# library which uses logging
#set(EASYLOGGINGPP_ENABLE $<TARGET_OBJECTS:easyloggingpp_object> CACHE STRING
#	"Link with this object to setup easyloggingpp automatically")

set(ELPP_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
