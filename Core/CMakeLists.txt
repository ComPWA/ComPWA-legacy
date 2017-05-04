# Create ComPWA Core library.
FILE(GLOB lib_srcs *.cpp)
FILE(GLOB lib_headers *.hpp)
#FILE(GLOB lib_ignore *Test.cpp)
list(REMOVE_ITEM lib_srcs "${CMAKE_CURRENT_SOURCE_DIR}/ParameterListTest.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/ParameterTest.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/Core_ext.cpp")

add_library ( ComPWA_Core SHARED ${lib_srcs} ${lib_headers} )

target_link_libraries( ComPWA_Core ${Boost_LIBRARIES} )

target_include_directories ( ComPWA_Core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

#install (TARGETS ComPWA_Core DESTINATION lib/ComPWA)
install (FILES ${lib_headers} DESTINATION include/ComPWA)

#find_package( SWIG )
#IF(SWIG_FOUND)
#    INCLUDE(${SWIG_USE_FILE})
#
#    find_package(PythonLibs)
#    IF(PYTHONLIBS_FOUND)
#        INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
#
#        INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
#
#        SET(CMAKE_SWIG_FLAGS "")
#
#        SET_SOURCE_FILES_PROPERTIES(CoreToScript.i PROPERTIES CPLUSPLUS ON)
#        #SET_SOURCE_FILES_PROPERTIES(CoreToScript.i PROPERTIES SWIG_FLAGS "-includeall")
#        SWIG_ADD_MODULE(core python CoreToScript.i)
#        SWIG_LINK_LIBRARIES(core ${PYTHON_LIBRARIES})
#    ENDIF()
#ENDIF()

IF(PYTHONLIBS_FOUND)
    # Define the wrapper library that wraps our library
    add_library( Core_ext SHARED Core_ext.cpp )
    target_link_libraries( Core_ext ${Boost_LIBRARIES} ComPWA_Core )
    # don't prepend wrapper library name with lib
    set_target_properties( Core_ext PROPERTIES PREFIX "" )
ENDIF()
