########################
# Create qft++ library #
########################

set(SOURCE_LIST src/TensorIndex.C src/Tensor.C)

set(HEADER_LIST
    ${CMAKE_INSTALL_INCLUDEDIR}/qft++/MatrixOfTensors.h
    ${CMAKE_INSTALL_INCLUDEDIR}/qft++/SpecialTensors.h
    ${CMAKE_INSTALL_INCLUDEDIR}/qft++/Tensor.h
    ${CMAKE_INSTALL_INCLUDEDIR}/qft++/Tensor_Base.h
    ${CMAKE_INSTALL_INCLUDEDIR}/qft++/TensorIndex.h
    ${CMAKE_INSTALL_INCLUDEDIR}/qft++/Vector4.h
    ${CMAKE_INSTALL_INCLUDEDIR}/qft++/WignerD.h
    ${CMAKE_INSTALL_INCLUDEDIR}/qft++/OperationType.h
    ${CMAKE_INSTALL_INCLUDEDIR}/qft++/Type.h
    ${CMAKE_INSTALL_INCLUDEDIR}/qft++/Matrix_Base.h
    ${CMAKE_INSTALL_INCLUDEDIR}/qft++/SelectiveInclusion.h
    ${CMAKE_INSTALL_INCLUDEDIR}/qft++/Conversion.h
)

add_library(qft++ ${SOURCE_LIST} ${HEADER_LIST})
target_include_directories(qft++ PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(QFTPP_INCLUDE_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    CACHE INTERNAL ""
)

#
# Install
#
install(FILES ${HEADER_LIST}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ThirdParty/qft++
)
install(
  TARGETS qft++
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
