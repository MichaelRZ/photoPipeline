// #define Vec3 Vec3_mvg
#include "mvg.h"
// #undef Vec3
// #undef Load
#include "mvs.h"
#include <iostream>
#include <vector>
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
// #undef PLANAR
// #include "colmap.h"

int main(int argc, char **argv)
{
  //OpenMVG SECTION ......................................................................

  // std::cout << OPENMVG_USE_OPENMP << std::endl;
  std::string sImageDir = "/Users/michaelr/Desktop/imitation/images/",
  sfileDatabase = "/Users/michaelr/Desktop/imitation/sensor_width_camera_database.txt",
  sOutputDir = "/Users/michaelr/Desktop/imitation/mvgOutputs/";

  //0. Intrinsics analysis             openMVG_main_SfMInit_ImageListing
  initImageList(sImageDir, sOutputDir, sfileDatabase);
  //1. Compute features                openMVG_main_ComputeFeatures  
  // getFeatures(sOutputDir+"sfm_data.json", sOutputDir);
  // //2. Compute pairs                   openMVG_main_PairGenerator
  // imagePairer(sOutputDir+"sfm_data.json", sOutputDir+"pairs.bin");
  // //3. Compute matches                 openMVG_main_ComputeMatches
  // //BREAKS UNDER OPENMP...
  // imageMatcher(sOutputDir+"sfm_data.json", sOutputDir+"pairMatches.bin", sOutputDir+"pairs.bin");
  // //4. Filter matches                  openMVG_main_GeometricFilter
  // imageFilterGeo(sOutputDir+"sfm_data.json", sOutputDir+"pairMatches.bin", sOutputDir+"matches.f.bin");
  // //5. Incremental reconstruction      openMVG_main_SfM
  // imageSfm(sOutputDir+"sfm_data.json", sOutputDir, sOutputDir+"sfm/");

  // if (!stlplus::folder_exists(sOutputDir+"mvs/"))
  // {
  //   if (!stlplus::folder_create(sOutputDir+"mvs/"))
  //   {
  //     std::cout << "Cannot create the output directory\n";
  //     return 1;
  //   }
  // }
  // //11. Export to openMVS              openMVG_main_openMVG2openMVS
  // mvgToMvs(sOutputDir+"sfm/sfm_data.bin", sOutputDir+"mvs/scene.mvs", sOutputDir+"mvs/images/");

  //COLMAP SECTION .....................................................................

  // // output_path: pathlib.Path
  // // image_dir: pathlib.Path

  // // output_path.mkdir()
  // // mvs_path = output_path / "mvs"
  // // database_path = output_path / "database.db"

  // // pycolmap.extract_features(database_path, image_dir)
  // ExtractFeatures(sImageDir+"features.db",sImageDir);
  // // pycolmap.match_exhaustive(database_path)
  // // maps = pycolmap.incremental_mapping(database_path, image_dir, output_path)
  // // maps[0].write(output_path)
  // // # dense reconstruction
  // // pycolmap.undistort_images(mvs_path, output_path, image_dir)
  // // pycolmap.patch_match_stereo(mvs_path)  # requires compilation with CUDA
  // // pycolmap.stereo_fusion(mvs_path / "dense.ply", mvs_path)

  // ["Feature Extractor",            # 12
  //   COLMAP_BIN,
  //   ["feature_extractor", "--database_path", "%matches_dir%"+FOLDER_DELIM+"database.db", "--image_path", "%input_dir%"]],
  // std::vector<std::string> arguments_COL= {"feature_extractor", "--database_path",sOutputDir+"database.db", 
  // "--image_path", sImageDir};
  // colmapCall(arguments_COL);
  // ["Exhaustive Matcher",           # 13
  //   COLMAP_BIN,
  //   ["exhaustive_matcher", "--database_path", "%matches_dir%"+FOLDER_DELIM+"database.db"]],
  // ["Mapper",                       # 14
  //   COLMAP_BIN,
  //   ["mapper", "--database_path", "%matches_dir%"+FOLDER_DELIM+"database.db", "--image_path", "%input_dir%", "--output_path", "%reconstruction_dir%"]],
  // ["Image Undistorter",            # 15
  //   COLMAP_BIN,
  //   ["image_undistorter", "--image_path", "%input_dir%", "--input_path", "%reconstruction_dir%"+FOLDER_DELIM+"0", "--output_path", "%reconstruction_dir%"+FOLDER_DELIM+"dense", "--output_type", "COLMAP"]],
  // ["Export to openMVS",            # 16
  //   os.path.join(OPENMVS_BIN, "InterfaceCOLMAP"),
  //   ["-i", "%reconstruction_dir%"+FOLDER_DELIM+"dense", "-o", "scene.mvs", "--image-folder", "%reconstruction_dir%"+FOLDER_DELIM+"dense"+FOLDER_DELIM+"images", "-w", "\"%mvs_dir%\""]],
            


  // 17. Densify point-cloud            DensifyPointCloud
  // std::vector<std::string> arguments = {"--dir", "/some_path"};
  // os.path.join(OPENMVS_BIN, "DensifyPointCloud"),
  // ["scene.mvs", "--dense-config-file", "Densify.ini", "--resolution-level", "1", "--number-views", "8", "-w", "\"%mvs_dir%\""]],
  std::vector<std::string> arguments = {"densify", "scene.mvs","--dense-config-file", 
  "Densify.ini", "--resolution-level", "1", "--number-views", "8", "-w", sOutputDir+"mvs/"};
  // densePtCld(arguments);

  // 18. Reconstruct the mesh           ReconstructMesh
  //Breaks from openmp, added #undef _USE_OPENMP to mvs.cpp to fix on 2021 version
  //might not?
  //Nevermind, only works in "Debug" mode.
  arguments = {"mesh" ,"scene_dense.mvs", "-o", "scene_dense.ply", "-w", sOutputDir+"mvs/"};
  reconstructMesh(arguments);

  // 19. Refine the mesh                RefineMesh


  // 20. Texture the mesh               TextureMesh
  //unneeded
  return 0;
}
