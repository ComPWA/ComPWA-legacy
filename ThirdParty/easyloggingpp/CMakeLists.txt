#
# Easyloggingpp - Logging library
#

# Create object file
add_library(easyloggingpp_obj OBJECT easylogging++.cc)

target_compile_options(easyloggingpp_obj PUBLIC -DELPP_NO_DEFAULT_LOG_FILE)

set(ELPP_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
install( FILES easylogging++.h DESTINATION include/)
