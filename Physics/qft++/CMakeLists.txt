# Create qft++ target.

SET( lib_srcs TensorIndex.C Tensor.C )
SET( lib_headers  MatrixOfTensors.h SpecialTensors.h Tensor.h Tensor_Base.h
      TensorIndex.h Vector4.h WignerD.h OperationType.h Type.h Matrix_Base.h
      SelectiveInclusion.h Conversion.h )

add_library ( qft++ SHARED ${lib_srcs} ${lib_headers})
#target_link_libraries( qft++ ${Boost_LIBRARIES} )
#target_include_directories( HelicityFormalism PUBLIC ${Boost_INCLUDE_DIR} )
install (FILES ${lib_headers} DESTINATION include/ComPWA)

