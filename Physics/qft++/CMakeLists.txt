# Create qft++ target.

SET( lib_srcs TensorIndex.C )
SET( lib_headers  MatrixOfTensors.h SpecialTensors.h Tensor.h Tensor_Base.h
      TensorIndex.h Vector4.h WignerD.h OperationType.h Type.h Matrix_Base.h
      SelectiveInclusion.h )

add_library ( ComPWA_qft++ SHARED ${lib_srcs} ${lib_headers})
#target_link_libraries( ComPWA_qft++ ${Boost_LIBRARIES} )
#target_include_directories( ComPWA_HelicityFormalism PUBLIC ${Boost_INCLUDE_DIR} )
install (FILES ${lib_headers} DESTINATION include/ComPWA)

