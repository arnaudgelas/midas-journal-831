set(PROJECT_TEST_PATH ${EXECUTABLE_OUTPUT_PATH})
set(TEST_DATA_ROOT ${CMAKE_SOURCE_DIR}/Data)
set(TEST_OUTPUT ${CMAKE_BINARY_DIR}/Testing/Temporary)

if("${ITK_VERSION_MAJOR}" LESS 4)
  set(USED_LIBS ITKIO ITKNumerics ITKQuadEdgeMesh)
else("${ITK_VERSION_MAJOR}" LESS 4)
  set(USED_LIBS ${ITK_LIBRARIES})
endif("${ITK_VERSION_MAJOR}" LESS 4)

add_executable(SubdivisionQuadEdgeMeshFilterTest SubdivisionQuadEdgeMeshFilterTest.cxx )
target_link_libraries(SubdivisionQuadEdgeMeshFilterTest ${USED_LIBS})

add_test(LinearTriangleCellSubdivisionQuadEdgeMeshFilterTest_1
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/ico.vtk
  ${TEST_OUTPUT}/ico_linear.vtk
  1
  1
  )

add_test(LinearTriangleCellSubdivisionQuadEdgeMeshFilterTest_2
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/ico.vtk
  ${TEST_OUTPUT}/ico_linear_2.vtk
  1
  2
  )

add_test(ModifiedButterflyTriangleCellSubdivisionQuadEdgeMeshFilterTest_1
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/ico.vtk
  ${TEST_OUTPUT}/ico_butterfly.vtk
  0
  1
  )

add_test(ModifiedButterflyTriangleCellSubdivisionQuadEdgeMeshFilterTest_2
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/ico.vtk
  ${TEST_OUTPUT}/ico_butterfly_2.vtk
  0
  2
  )

add_test(LoopTriangleCellSubdivisionQuadEdgeMeshFilterTest_1
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/ico.vtk
  ${TEST_OUTPUT}/ico_loop.vtk
  2
  1
  )

add_test(LoopTriangleCellSubdivisionQuadEdgeMeshFilterTest_2
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/ico.vtk
  ${TEST_OUTPUT}/ico_loop_2.vtk
  2
  2
  )

add_test(SquareThreeTriangleCellSubdivisionQuadEdgeMeshFilterTest_1
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/ico.vtk
  ${TEST_OUTPUT}/ico_square3.vtk
  3
  1
  )

add_test(SquareThreeTriangleCellSubdivisionQuadEdgeMeshFilterTest_2
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/ico.vtk
  ${TEST_OUTPUT}/ico_square3_2.vtk
  3
  2
  )

add_test(AdaptiveLinearTriangleCellSubdivisionQuadEdgeMeshFilterTest_1
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/ico.vtk
  ${TEST_OUTPUT}/ico_alinear.vtk
  1
  1
  1
  )

add_test(AdaptiveLinearTriangleCellSubdivisionQuadEdgeMeshFilterTest_2
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/ico.vtk
  ${TEST_OUTPUT}/ico_alinear_2.vtk
  1
  2
  1
  )

add_test(AdaptiveModifiedButterflyTriangleCellSubdivisionQuadEdgeMeshFilterTest_1
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/ico.vtk
  ${TEST_OUTPUT}/ico_abutterfly.vtk
  0
  1
  1
  )

add_test(AdaptiveModifiedButterflyTriangleCellSubdivisionQuadEdgeMeshFilterTest_2
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/ico.vtk
  ${TEST_OUTPUT}/ico_abutterfly_2.vtk
  0
  2
  1
  )

add_test(AdaptiveLoopTriangleCellSubdivisionQuadEdgeMeshFilterTest_1
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/ico.vtk
  ${TEST_OUTPUT}/ico_aloop.vtk
  2
  1
  1
  )

add_test(AdaptiveLoopTriangleCellSubdivisionQuadEdgeMeshFilterTest_2
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/ico.vtk
  ${TEST_OUTPUT}/ico_aloop_2.vtk
  2
  2
  1
  )

add_test(AdaptiveSquareThreeTriangleCellSubdivisionQuadEdgeMeshFilterTest_1
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/ico.vtk
  ${TEST_OUTPUT}/ico_asquare3.vtk
  3
  1
  1
  )

add_test(AdaptiveSquareThreeTriangleCellSubdivisionQuadEdgeMeshFilterTest_2
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/ico.vtk
  ${TEST_OUTPUT}/ico_asquare3_2.vtk
  3
  2
  1
  )

add_test(LinearTriangleCellSubdivisionQuadEdgeMeshFilterPlaneTest_1
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/plane.vtk
  ${TEST_OUTPUT}/plane_linear.vtk
  1
  1
  )

add_test(LinearTriangleCellSubdivisionQuadEdgeMeshFilterPlaneTest_2
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/plane.vtk
  ${TEST_OUTPUT}/plane_linear_2.vtk
  1
  2
  )

add_test(ModifiedButterflyTriangleCellSubdivisionQuadEdgeMeshFilterPlaneTest_1
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/plane.vtk
  ${TEST_OUTPUT}/plane_butterfly.vtk
  0
  1
  )

add_test(ModifiedButterflyTriangleCellSubdivisionQuadEdgeMeshFilterPlaneTest_2
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/plane.vtk
  ${TEST_OUTPUT}/plane_butterfly_2.vtk
  0
  2
  )

add_test(LoopTriangleCellSubdivisionQuadEdgeMeshFilterPlaneTest_1
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/plane.vtk
  ${TEST_OUTPUT}/plane_loop.vtk
  2
  1
  )

add_test(LoopTriangleCellSubdivisionQuadEdgeMeshFilterPlaneTest_2
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/plane.vtk
  ${TEST_OUTPUT}/plane_loop_2.vtk
  2
  2
  )

add_test(SquareThreeTriangleCellSubdivisionQuadEdgeMeshFilterPlaneTest_1
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/plane.vtk
  ${TEST_OUTPUT}/plane_square3.vtk
  3
  1
  )

add_test(SquareThreeTriangleCellSubdivisionQuadEdgeMeshFilterPlaneTest_2
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/plane.vtk
  ${TEST_OUTPUT}/plane_square3_2.vtk
  3
  2
  )

add_test(AdaptiveLinearTriangleCellSubdivisionQuadEdgeMeshFilterPlaneTest_1
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/plane.vtk
  ${TEST_OUTPUT}/plane_alinear.vtk
  1
  1
  1
  )

add_test(AdaptiveLinearTriangleCellSubdivisionQuadEdgeMeshFilterPlaneTest_2
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/plane.vtk
  ${TEST_OUTPUT}/plane_alinear_2.vtk
  1
  2
  1
  )

add_test(AdaptiveModifiedButterflyTriangleCellSubdivisionQuadEdgeMeshFilterPlaneTest_1
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/plane.vtk
  ${TEST_OUTPUT}/plane_abutterfly.vtk
  0
  1
  1
  )

add_test(AdaptiveModifiedButterflyTriangleCellSubdivisionQuadEdgeMeshFilterPlaneTest_2
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/plane.vtk
  ${TEST_OUTPUT}/plane_abutterfly_2.vtk
  0
  2
  1
  )

add_test(AdaptiveLoopTriangleCellSubdivisionQuadEdgeMeshFilterPlaneTest_1
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/plane.vtk
  ${TEST_OUTPUT}/plane_aloop.vtk
  2
  1
  0
  )

add_test(AdaptiveLoopTriangleCellSubdivisionQuadEdgeMeshFilterPlaneTest_2
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/plane.vtk
  ${TEST_OUTPUT}/plane_aloop_2.vtk
  2
  2
  0
  )

add_test(AdaptiveSquareThreeTriangleCellSubdivisionQuadEdgeMeshFilterPlaneTest_1
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/plane.vtk
  ${TEST_OUTPUT}/plane_asquare3.vtk
  3
  1
  0
  )

add_test(AdaptiveSquareThreeTriangleCellSubdivisionQuadEdgeMeshFilterPlaneTest_2
  ${PROJECT_TEST_PATH}/SubdivisionQuadEdgeMeshFilterTest
  ${TEST_DATA_ROOT}/plane.vtk
  ${TEST_OUTPUT}/plane_asquare3_2.vtk
  3
  2
  0
  )
