#
# Easyloggingpp - Logging library
#

# Create object file
add_library(easyloggingpp_obj OBJECT easylogging++.cc)

set(ELPP_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")

install (FILES easylogging++.h
    DESTINATION include/ThirdParty/easyloggingpp/
)
