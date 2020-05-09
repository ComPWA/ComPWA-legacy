########################
# Create qft++ library #
########################

set(SOURCE_LIST src/TensorIndex.C src/Tensor.C)

set(HEADER_LIST
    include/qft++/MatrixOfTensors.h
    include/qft++/SpecialTensors.h
    include/qft++/Tensor.h
    include/qft++/Tensor_Base.h
    include/qft++/TensorIndex.h
    include/qft++/Vector4.h
    include/qft++/WignerD.h
    include/qft++/OperationType.h
    include/qft++/Type.h
    include/qft++/Matrix_Base.h
    include/qft++/SelectiveInclusion.h
    include/qft++/Conversion.h)

add_library(qft++ ${SOURCE_LIST} ${HEADER_LIST})
target_include_directories(qft++ PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(QFTPP_INCLUDE_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    CACHE INTERNAL "")

#
# Install
#
install(FILES ${HEADER_LIST} DESTINATION include/ThirdParty/qft++)
install(
  TARGETS qft++
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
