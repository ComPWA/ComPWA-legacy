########################
# Create qft++ library #
########################

SET( lib_srcs
  src/TensorIndex.C
  src/Tensor.C
)

SET( lib_headers
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
  include/qft++/Conversion.h
)

add_library( qft++
  OBJECT ${lib_srcs} ${lib_headers}
)
target_include_directories( qft++ PUBLIC 
  ${CMAKE_CURRENT_SOURCE_DIR}/include )

set(QFTPP $<TARGET_OBJECTS:qft++> CACHE STRING
  "Link with this object to qft++")

set(QFTPP_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include CACHE INTERNAL "")
#
# Install
#
install (FILES ${lib_headers}
    DESTINATION include/ComPWA/qft++
)

#install(TARGETS qft++
#        LIBRARY DESTINATION lib/ComPWA
#)
