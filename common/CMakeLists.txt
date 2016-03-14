set(srcs
  Data_Obj/Image3D.cpp
  Data_Obj/TriangleMesh.cpp
  Data_Obj/Triplet.cpp

  IO/LoadBigImage.cpp
  IO/LoadImage3DMPI.cpp

  Reporting/Timer.cpp
  Reporting/YAML_Doc.cpp
  Reporting/YAML_Element.cpp

  Constants/MarchingCubesTables.cpp
  Constants/TypeInfo.cpp

  User_Interface/UI.cpp

  Algorithm/BlockMarchFunctor.cpp
  Algorithm/EdgeIndexer.cpp
  Algorithm/gradients.cpp
  Algorithm/Image3DReader.cpp
  Algorithm/MarchAlgorithm.cpp

  GeneralContext/GeneralContext.cpp
  )

add_library(mantevo_marching_cubes_common
  ${srcs}
  )
