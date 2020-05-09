set(PACKAGE_NAME EvtGen)

# Define library
set(SOURCE_LIST
    EvtBlattWeisskopf.cpp
    EvtComplex.cpp
    EvtConst.cpp
    EvtCyclic3.cpp
    EvtDalitzCoord.cpp
    EvtDalitzPlot.cpp
    EvtDalitzPoint.cpp
    EvtDalitzReso.cpp
    EvtFlatte.cpp
    EvtGenKine.cpp
    EvtRandom.cpp
    EvtReport.cpp
    EvtSpinType.cpp
    EvtTensor3C.cpp
    EvtTensor4C.cpp
    EvtTwoBodyKine.cpp
    EvtTwoBodyVertex.cpp
    EvtVector3C.cpp
    EvtVector3R.cpp
    EvtVector4C.cpp
    EvtVector4R.cpp
    EvtdFunction.cpp
    EvtdFunctionSingle.cpp
)
add_library(${PACKAGE_NAME} ${SOURCE_LIST})
target_include_directories(
  ${PACKAGE_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
                         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Install and export
set(HEADER_LIST
    EvtBlattWeisskopf.hh
    EvtComplex.hh
    EvtConst.hh
    EvtCyclic3.hh
    EvtDalitzCoord.hh
    EvtDalitzPlot.hh
    EvtDalitzPoint.hh
    EvtDalitzReso.hh
    EvtFlatte.hh
    EvtGenKine.hh
    EvtMacros.hh
    EvtMatrix.hh
    EvtPatches.hh
    EvtRandom.hh
    EvtRandomEngine.hh
    EvtReport.hh
    EvtSpinType.hh
    EvtTensor3C.hh
    EvtTensor4C.hh
    EvtTwoBodyKine.hh
    EvtTwoBodyVertex.hh
    EvtVector3C.hh
    EvtVector3R.hh
    EvtVector4C.hh
    EvtVector4R.hh
    EvtdFunction.hh
    EvtdFunctionSingle.hh
)
install(FILES ${HEADER_LIST}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PACKAGE_NAME}
)
install(
  TARGETS ${PACKAGE_NAME}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
