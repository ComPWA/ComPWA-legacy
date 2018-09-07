add_subdirectory(Minuit2)
add_subdirectory(Geneva)

install (FILES Optimizer.hpp
    DESTINATION include/Optimizer
)

#
# TESTING
#
# tests are not implemented correctly atm
return()

if(${ROOT_MINUIT2_FOUND})
add_executable(minuitTestApp test/MinuitTestApp.cpp)
target_link_libraries(minuitTestApp
	Core
    ROOT::Minuit2
)
set_target_properties(minuitTestApp
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
)
add_test(NAME minuitTestApp
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/minuitTestApp
)
elseif(${GENEVA_FOUND})
add_executable(genevaTestApp test/GenevaTestApp.cpp)
target_link_libraries(genevaTestApp
	Core
    ${GENEVA_LIBRARIES}
)
set_target_properties(genevaTestApp
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
)
add_test(NAME genevaTestApp
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/genevaTestApp
)
elseif(${ROOT_MINUIT2_FOUND} AND ${GENEVA_FOUND})
add_executable(optimizerTestApp test/OptimizerTestApp.cpp)
target_link_libraries(optimizerTestApp
	Core
    ROOT::Minuit2
    ${GENEVA_LIBRARIES}
)
set_target_properties(optimizerTestApp
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
)
add_test(NAME optimizerTestApp
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/bin/test/
    COMMAND ${PROJECT_BINARY_DIR}/bin/test/optimizerTestApp
)
endif()
