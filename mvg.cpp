#include "mvg.h"

// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013, 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//initlist requirements
#include "openMVG/cameras/cameras.hpp"
#include "openMVG/exif/exif_IO_EasyExif.hpp"
#include "openMVG/exif/sensor_width_database/ParseDatabase.hpp"
#include "openMVG/geodesy/geodesy.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_utils.hpp"
#include "openMVG/sfm/sfm_view.hpp"
#include "openMVG/sfm/sfm_view_priors.hpp"
#include "openMVG/system/loggerprogress.hpp"
#include "openMVG/types.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <fstream>
#include <memory>
#include <string>
#include <utility>

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::exif;
using namespace openMVG::geodesy;
using namespace openMVG::image;
using namespace openMVG::sfm;

//feature get requirements
// The <cereal/archives> headers are special and must be included first.
#include <cereal/archives/json.hpp>

#include "openMVG/features/akaze/image_describer_akaze_io.hpp"

#include "openMVG/features/sift/SIFT_Anatomy_Image_Describer_io.hpp"
#include "openMVG/features/regions_factory_io.hpp"
#include "openMVG/system/logger.hpp"
#include "openMVG/system/timer.hpp"

#include "nonFree/sift/SIFT_describer_io.hpp"

#include <cereal/details/helpers.hpp>

#include <atomic>
#include <cstdlib>
#include <fstream>

#ifdef OPENMVG_USE_OPENMP
#include <omp.h>
#endif

using namespace openMVG::features;


//pair gen requirements
#include "openMVG/matching_image_collection/Pair_Builder.hpp"

#include <iostream>

/**
 * @brief Current list of available pair mode
 *
 */
enum EPairMode
{
  PAIR_EXHAUSTIVE = 0, // Build every combination of image pairs
  PAIR_CONTIGUOUS = 1  // Only consecutive image pairs (useful for video mode)
};


//compute matches requirements
#include "openMVG/graph/graph.hpp"
#include "openMVG/graph/graph_stats.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/matching/pairwiseAdjacencyDisplay.hpp"
#include "openMVG/matching_image_collection/Cascade_Hashing_Matcher_Regions.hpp"
#include "openMVG/matching_image_collection/Matcher_Regions.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_preemptive_regions_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider_cache.hpp"
#include "openMVG/stl/stl.hpp"

using namespace openMVG::matching;
using namespace openMVG::matching_image_collection;

//geometric filter requirements
#include "openMVG/features/akaze/image_describer_akaze.hpp"
#include "openMVG/features/descriptor.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/matching_image_collection/E_ACRobust.hpp"
#include "openMVG/matching_image_collection/E_ACRobust_Angular.hpp"
#include "openMVG/matching_image_collection/Eo_Robust.hpp"
#include "openMVG/matching_image_collection/F_ACRobust.hpp"
#include "openMVG/matching_image_collection/GeometricFilter.hpp"
#include "openMVG/matching_image_collection/H_ACRobust.hpp"

#include <locale>

using namespace openMVG::robust;

//sfm requirements
#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/cameras/Cameras_Common_command_line_helper.hpp"

#include "openMVG/sfm/pipelines/sfm_matches_provider.hpp"
#include "openMVG/sfm/sfm_report.hpp"

// SfM Engines
#include "openMVG/sfm/pipelines/sequential/sequential_SfM.hpp"
#include "openMVG/sfm/pipelines/sequential/sequential_SfM2.hpp"
#include "openMVG/sfm/pipelines/sequential/SfmSceneInitializerMaxPair.hpp"
#include "openMVG/sfm/pipelines/sequential/SfmSceneInitializerStellar.hpp"
#include "openMVG/sfm/pipelines/global/GlobalSfM_rotation_averaging.hpp"
#include "openMVG/sfm/pipelines/global/GlobalSfM_translation_averaging.hpp"
#include "openMVG/sfm/pipelines/global/sfm_global_engine_relative_motions.hpp"


//mvg to mvs requirements
#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/cameras/Camera_undistort_image.hpp"

// #define _USE_EIGEN
#include "InterfaceMVS.h"

using namespace openMVG::geometry;


/// Check that Kmatrix is a string like "f;0;ppx;0;f;ppy;0;0;1"
/// With f,ppx,ppy as valid numerical value
bool checkIntrinsicStringValidity(const std::string & Kmatrix, double & focal, double & ppx, double & ppy)
{
  std::vector<std::string> vec_str;
  stl::split(Kmatrix, ';', vec_str);
  if (vec_str.size() != 9)  {
    OPENMVG_LOG_ERROR << "\n Missing ';' character";
    return false;
  }
  // Check that all K matrix value are valid numbers
  for (size_t i = 0; i < vec_str.size(); ++i) {
    double readvalue = 0.0;
    std::stringstream ss;
    ss.str(vec_str[i]);
    if (! (ss >> readvalue) )  {
      OPENMVG_LOG_ERROR << "\n Used an invalid not a number character";
      return false;
    }
    if (i==0) focal = readvalue;
    if (i==2) ppx = readvalue;
    if (i==5) ppy = readvalue;
  }
  return true;
}

bool getGPS
(
  const std::string & filename,
  const int & GPS_to_XYZ_method,
  Vec3 & pose_center
)
{
  std::unique_ptr<Exif_IO> exifReader(new Exif_IO_EasyExif);
  if (exifReader)
  {
    // Try to parse EXIF metada & check existence of EXIF data
    if ( exifReader->open( filename ) && exifReader->doesHaveExifInfo() )
    {
      // Check existence of GPS coordinates
      double latitude, longitude, altitude;
      if ( exifReader->GPSLatitude( &latitude ) &&
           exifReader->GPSLongitude( &longitude ) &&
           exifReader->GPSAltitude( &altitude ) )
      {
        // Add ECEF or UTM XYZ position to the GPS position array
        switch (GPS_to_XYZ_method)
        {
          case 1:
            pose_center = lla_to_utm( latitude, longitude, altitude );
            break;
          case 0:
          default:
            pose_center = lla_to_ecef( latitude, longitude, altitude );
            break;
        }
        return true;
      }
    }
  }
  return false;
}


/// Check string of prior weights
std::pair<bool, Vec3> checkPriorWeightsString
(
  const std::string &sWeights
)
{
  std::pair<bool, Vec3> val(true, Vec3::Zero());
  std::vector<std::string> vec_str;
  stl::split(sWeights, ';', vec_str);
  if (vec_str.size() != 3)
  {
    OPENMVG_LOG_ERROR << "Missing ';' character in prior weights";
    val.first = false;
  }
  // Check that all weight values are valid numbers
  for (size_t i = 0; i < vec_str.size(); ++i)
  {
    double readvalue = 0.0;
    std::stringstream ss;
    ss.str(vec_str[i]);
    if (! (ss >> readvalue) )  {
      OPENMVG_LOG_ERROR << "Used an invalid not a number character in local frame origin";
      val.first = false;
    }
    val.second[i] = readvalue;
  }
  return val;
}
void GroupSharedIntrinsics2(SfM_Data & sfm_data)
{
  Views & views = sfm_data.views;
  Intrinsics & intrinsics = sfm_data.intrinsics;

  // Build hash & build a set of the hash in order to maintain unique Ids
  std::set<size_t> hash_index;
  std::vector<size_t> hash_value;

  for (const auto & intrinsic_it : intrinsics)
  {
    const cameras::IntrinsicBase * intrinsicData = intrinsic_it.second.get();
    const size_t hashVal = intrinsicData->hashValue();
    hash_index.insert(hashVal);
    hash_value.push_back(hashVal);
  }

  // From hash_value(s) compute the new index (old to new indexing)
  Hash_Map<IndexT, IndexT> old_new_reindex;
  size_t i = 0;
  for (const auto & intrinsic_it : intrinsics)
  {
    old_new_reindex[intrinsic_it.first] = std::distance(hash_index.cbegin(), hash_index.find(hash_value[i]));
    ++i;
  }
  //--> Save only the required Intrinsics (do not need to keep all the copy)
  Intrinsics intrinsic_updated;
  for (const auto & intrinsic_it : intrinsics)
  {
    intrinsic_updated[old_new_reindex[intrinsic_it.first]] = intrinsics[intrinsic_it.first];
  }
  // Update intrinsics (keep only the necessary ones) -> swapping
  intrinsics.swap(intrinsic_updated);

  // Update views intrinsic IDs (since some intrinsic position have changed in the map)
  for (auto & view_it: views)
  {
    View * v = view_it.second.get();
    // Update the Id only if a corresponding index exists
    if (old_new_reindex.count(v->id_intrinsic))
      v->id_intrinsic = old_new_reindex[v->id_intrinsic];
  }
}

//
// Create the description of an input image dataset for OpenMVG toolsuite
// - Export a SfM_Data file with View & Intrinsic data
//
int initImageList(std::string sImageDir, std::string sOutputDir, std::string sfileDatabase){
  std::string sKmatrix;

  std::string sPriorWeights = "1.0;1.0;1.0";
  std::pair<bool, Vec3> prior_w_info(false, Vec3());

  int i_User_camera_model = PINHOLE_CAMERA_RADIAL3;

  bool b_Group_camera_model = true;

  int i_GPS_XYZ_method = 0;

  double focal_pixels = -1.0;

  const bool b_Use_pose_prior = false;//cmd.used('P');
  OPENMVG_LOG_INFO << " You called : " << " "//argv[0]
    << "\n--imageDirectory " << sImageDir
    << "\n--sensorWidthDatabase " << sfileDatabase
    << "\n--outputDirectory " << sOutputDir
    << "\n--focal " << focal_pixels
    << "\n--intrinsics " << sKmatrix
    << "\n--camera_model " << i_User_camera_model
    << "\n--group_camera_model " << b_Group_camera_model
    << "\n--use_pose_prior " << b_Use_pose_prior
    << "\n--prior_weights " << sPriorWeights
    << "\n--gps_to_xyz_method " << i_GPS_XYZ_method;

  // Expected properties for each image
  double width = -1, height = -1, focal = -1, ppx = -1,  ppy = -1;

  const EINTRINSIC e_User_camera_model = EINTRINSIC(i_User_camera_model);

  if ( !stlplus::folder_exists( sImageDir ) )
  {
    OPENMVG_LOG_ERROR << "The input directory doesn't exist";
    return EXIT_FAILURE;
  }

  if (sOutputDir.empty())
  {
    OPENMVG_LOG_ERROR << "Invalid output directory";
    return EXIT_FAILURE;
  }

  if ( !stlplus::folder_exists( sOutputDir ) )
  {
    if ( !stlplus::folder_create( sOutputDir ))
    {
      OPENMVG_LOG_ERROR << "Cannot create output directory";
      return EXIT_FAILURE;
    }
  }

  if (sKmatrix.size() > 0 &&
    !checkIntrinsicStringValidity(sKmatrix, focal, ppx, ppy) )
  {
    OPENMVG_LOG_ERROR << "Invalid K matrix input";
    return EXIT_FAILURE;
  }

  if (sKmatrix.size() > 0 && focal_pixels != -1.0)
  {
    OPENMVG_LOG_ERROR << "Cannot combine -f and -k options";
    return EXIT_FAILURE;
  }

  std::vector<Datasheet> vec_database;
  if (!sfileDatabase.empty())
  {
    if ( !parseDatabase( sfileDatabase, vec_database ) )
    {
      OPENMVG_LOG_ERROR
       << "Invalid input database: " << sfileDatabase
       << ", please specify a valid file.";
      return EXIT_FAILURE;
    }
  }

  // Check if prior weights are given
  if (b_Use_pose_prior)
  {
    prior_w_info = checkPriorWeightsString(sPriorWeights);
  }

  std::vector<std::string> vec_image = stlplus::folder_files( sImageDir );
  std::sort(vec_image.begin(), vec_image.end());

  // Configure an empty scene with Views and their corresponding cameras
  SfM_Data sfm_data;
  sfm_data.s_root_path = sImageDir; // Setup main image root_path
  Views & views = sfm_data.views;
  Intrinsics & intrinsics = sfm_data.intrinsics;

  system::LoggerProgress my_progress_bar(vec_image.size(), "- Listing images -" );
  std::ostringstream error_report_stream;
  for ( std::vector<std::string>::const_iterator iter_image = vec_image.begin();
    iter_image != vec_image.end();
    ++iter_image, ++my_progress_bar )
  {
    // Read meta data to fill camera parameter (w,h,focal,ppx,ppy) fields.
    width = height = ppx = ppy = focal = -1.0;

    const std::string sImageFilename = stlplus::create_filespec( sImageDir, *iter_image );
    const std::string sImFilenamePart = stlplus::filename_part(sImageFilename);

    // Test if the image format is supported:
    if (openMVG::image::GetFormat(sImageFilename.c_str()) == openMVG::image::Unknown)
    {
      error_report_stream
          << sImFilenamePart << ": Unkown image file format." << "\n";
      continue; // image cannot be opened
    }

    if (sImFilenamePart.find("mask.png") != std::string::npos
       || sImFilenamePart.find("_mask.png") != std::string::npos)
    {
      error_report_stream
          << sImFilenamePart << " is a mask image" << "\n";
      continue;
    }

    ImageHeader imgHeader;
    if (!openMVG::image::ReadImageHeader(sImageFilename.c_str(), &imgHeader))
      continue; // image cannot be read

    width = imgHeader.width;
    height = imgHeader.height;
    ppx = width / 2.0;
    ppy = height / 2.0;


    // Consider the case where the focal is provided manually
    if (sKmatrix.size() > 0) // Known user calibration K matrix
    {
      if (!checkIntrinsicStringValidity(sKmatrix, focal, ppx, ppy))
        focal = -1.0;
    }
    else // User provided focal length value
      if (focal_pixels != -1 )
        focal = focal_pixels;

    // If not manually provided or wrongly provided
    if (focal == -1)
    {
      std::unique_ptr<Exif_IO> exifReader(new Exif_IO_EasyExif);
      exifReader->open( sImageFilename );

      const bool bHaveValidExifMetadata =
        exifReader->doesHaveExifInfo()
        && !exifReader->getModel().empty()
        && !exifReader->getBrand().empty();

      if (bHaveValidExifMetadata) // If image contains meta data
      {
        // Handle case where focal length is equal to 0
        if (exifReader->getFocal() == 0.0f)
        {
          error_report_stream
            << stlplus::basename_part(sImageFilename) << ": Focal length is missing." << "\n";
          focal = -1.0;
        }
        else
        // Create the image entry in the list file
        {
          const std::string sCamModel = exifReader->getBrand() + " " + exifReader->getModel();

          Datasheet datasheet;
          if ( getInfo( sCamModel, vec_database, datasheet ))
          {
            // The camera model was found in the database so we can compute it's approximated focal length
            const double ccdw = datasheet.sensorSize_;
            focal = std::max ( width, height ) * exifReader->getFocal() / ccdw;
          }
          else
          {
            error_report_stream
              << stlplus::basename_part(sImageFilename)
              << "\" model \"" << sCamModel << "\" doesn't exist in the database" << "\n"
              << "Please consider add your camera model and sensor width in the database." << "\n";
          }
        }
      }
    }
    // Build intrinsic parameter related to the view
    std::shared_ptr<IntrinsicBase> intrinsic;

    if (focal > 0 && ppx > 0 && ppy > 0 && width > 0 && height > 0)
    {
      // Create the desired camera type
      switch (e_User_camera_model)
      {
        case PINHOLE_CAMERA:
          intrinsic = std::make_shared<Pinhole_Intrinsic>
            (width, height, focal, ppx, ppy);
        break;
        case PINHOLE_CAMERA_RADIAL1:
          intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K1>
            (width, height, focal, ppx, ppy, 0.0); // setup no distortion as initial guess
        break;
        case PINHOLE_CAMERA_RADIAL3:
          intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K3>
            (width, height, focal, ppx, ppy, 0.0, 0.0, 0.0);  // setup no distortion as initial guess
        break;
        case PINHOLE_CAMERA_BROWN:
          intrinsic = std::make_shared<Pinhole_Intrinsic_Brown_T2>
            (width, height, focal, ppx, ppy, 0.0, 0.0, 0.0, 0.0, 0.0); // setup no distortion as initial guess
        break;
        case PINHOLE_CAMERA_FISHEYE:
          intrinsic = std::make_shared<Pinhole_Intrinsic_Fisheye>
            (width, height, focal, ppx, ppy, 0.0, 0.0, 0.0, 0.0); // setup no distortion as initial guess
        break;
        case CAMERA_SPHERICAL:
           intrinsic = std::make_shared<Intrinsic_Spherical>
             (width, height);
        break;
        default:
          OPENMVG_LOG_ERROR << "Error: unknown camera model: " << (int) e_User_camera_model;
          return EXIT_FAILURE;
      }
    }

    // Build the view corresponding to the image
    Vec3 pose_center;
    if (getGPS(sImageFilename, i_GPS_XYZ_method, pose_center) && b_Use_pose_prior)
    {
      ViewPriors v(*iter_image, views.size(), views.size(), views.size(), width, height);

      // Add intrinsic related to the image (if any)
      if (!intrinsic)
      {
        //Since the view have invalid intrinsic data
        // (export the view, with an invalid intrinsic field value)
        v.id_intrinsic = UndefinedIndexT;
      }
      else
      {
        // Add the defined intrinsic to the sfm_container
        intrinsics[v.id_intrinsic] = intrinsic;
      }

      v.b_use_pose_center_ = true;
      v.pose_center_ = pose_center;
      // prior weights
      if (prior_w_info.first == true)
      {
        v.center_weight_ = prior_w_info.second;
      }

      // Add the view to the sfm_container
      views[v.id_view] = std::make_shared<ViewPriors>(v);
    }
    else
    {
      View v(*iter_image, views.size(), views.size(), views.size(), width, height);

      // Add intrinsic related to the image (if any)
      if (!intrinsic)
      {
        //Since the view have invalid intrinsic data
        // (export the view, with an invalid intrinsic field value)
        v.id_intrinsic = UndefinedIndexT;
      }
      else
      {
        // Add the defined intrinsic to the sfm_container
        intrinsics[v.id_intrinsic] = intrinsic;
      }

      // Add the view to the sfm_container
      views[v.id_view] = std::make_shared<View>(v);
    }
  }

  // Display saved warning & error messages if any.
  if (!error_report_stream.str().empty())
  {
    OPENMVG_LOG_WARNING
      << "Warning & Error messages:\n"
      << error_report_stream.str();
  }
  std::cout << "group intrinsic" << std::endl;
  // Group camera that share common properties if desired (leads to more faster & stable BA).
  if (b_Group_camera_model)
  {
    GroupSharedIntrinsics2(sfm_data);
  }
  std::cout << "store stuff" << std::endl;
  // Store SfM_Data views & intrinsic data
  if (!Save(
    sfm_data,
    stlplus::create_filespec( sOutputDir, "sfm_data.json" ).c_str(),
    ESfM_Data(VIEWS|INTRINSICS)))
  {
    std::cout << "store failed" << std::endl;
    return EXIT_FAILURE;
  }

  OPENMVG_LOG_INFO
    << "SfMInit_ImageListing report:\n"
    << "listed #File(s): " << vec_image.size() << "\n"
    << "usable #File(s) listed in sfm_data: " << sfm_data.GetViews().size() << "\n"
    << "usable #Intrinsic(s) listed in sfm_data: " << sfm_data.GetIntrinsics().size();

  return EXIT_SUCCESS;
}

features::EDESCRIBER_PRESET stringToEnum(const std::string & sPreset)
{
  features::EDESCRIBER_PRESET preset;
  if (sPreset == "NORMAL")
    preset = features::NORMAL_PRESET;
  else
  if (sPreset == "HIGH")
    preset = features::HIGH_PRESET;
  else
  if (sPreset == "ULTRA")
    preset = features::ULTRA_PRESET;
  else
    preset = features::EDESCRIBER_PRESET(-1);
  return preset;
}
/// - Compute view image description (feature & descriptor extraction)
/// - Export computed data
int getFeatures(std::string sSfM_Data_Filename, std::string sOutDir)
{
  bool bUpRight = false;
  std::string sImage_Describer_Method = "SIFT";
  bool bForce = false;
  std::string sFeaturePreset = "";
  #ifdef OPENMVG_USE_OPENMP
    int iNumThreads = 0;
  #endif

  OPENMVG_LOG_INFO
    << " You called : " << "\n"
    //<< argv[0] << "\n"
    << "--input_file " << sSfM_Data_Filename << "\n"
    << "--outdir " << sOutDir << "\n"
    << "--describerMethod " << sImage_Describer_Method << "\n"
    << "--upright " << bUpRight << "\n"
    << "--describerPreset " << (sFeaturePreset.empty() ? "NORMAL" : sFeaturePreset) << "\n"
    << "--force " << bForce << "\n"
  #ifdef OPENMVG_USE_OPENMP
      << "--numThreads " << iNumThreads << "\n"
  #endif
    ;


  if (sOutDir.empty())
  {
    OPENMVG_LOG_ERROR << "\nIt is an invalid output directory";
    return EXIT_FAILURE;
  }

  // Create output dir
  if (!stlplus::folder_exists(sOutDir))
  {
    if (!stlplus::folder_create(sOutDir))
    {
      OPENMVG_LOG_ERROR << "Cannot create output directory";
      return EXIT_FAILURE;
    }
  }

  //---------------------------------------
  // a. Load input scene
  //---------------------------------------
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS))) {
    OPENMVG_LOG_ERROR
      << "The input file \""<< sSfM_Data_Filename << "\" cannot be read";
    return EXIT_FAILURE;
  }

  // b. Init the image_describer
  // - retrieve the used one in case of pre-computed features
  // - else create the desired one

  using namespace openMVG::features;
  std::unique_ptr<Image_describer> image_describer;

  const std::string sImage_describer = stlplus::create_filespec(sOutDir, "image_describer", "json");
  if (!bForce && stlplus::is_file(sImage_describer))
  {
    // Dynamically load the image_describer from the file (will restore old used settings)
    std::ifstream stream(sImage_describer.c_str());
    if (!stream)
      return EXIT_FAILURE;

    try
    {
      cereal::JSONInputArchive archive(stream);
      archive(cereal::make_nvp("image_describer", image_describer));
    }
    catch (const cereal::Exception & e)
    {
      OPENMVG_LOG_ERROR << e.what() << '\n'
        << "Cannot dynamically allocate the Image_describer interface.";
      return EXIT_FAILURE;
    }
  }
  else
  {
    // Create the desired Image_describer method.
    // Don't use a factory, perform direct allocation
    if (sImage_Describer_Method == "SIFT")
    {
      image_describer.reset(new SIFT_Image_describer
        (SIFT_Image_describer::Params(), !bUpRight));
    }
    else
    if (sImage_Describer_Method == "SIFT_ANATOMY")
    {
      image_describer.reset(
        new SIFT_Anatomy_Image_describer(SIFT_Anatomy_Image_describer::Params()));
    }
    else
    if (sImage_Describer_Method == "AKAZE_FLOAT")
    {
      image_describer = AKAZE_Image_describer::create
        (AKAZE_Image_describer::Params(AKAZE::Params(), AKAZE_MSURF), !bUpRight);
    }
    else
    if (sImage_Describer_Method == "AKAZE_MLDB")
    {
      image_describer = AKAZE_Image_describer::create
        (AKAZE_Image_describer::Params(AKAZE::Params(), AKAZE_MLDB), !bUpRight);
    }
    if (!image_describer)
    {
      OPENMVG_LOG_ERROR << "Cannot create the designed Image_describer:"
        << sImage_Describer_Method << ".";
      return EXIT_FAILURE;
    }
    else
    {
      if (!sFeaturePreset.empty())
      if (!image_describer->Set_configuration_preset(stringToEnum(sFeaturePreset)))
      {
        OPENMVG_LOG_ERROR << "Preset configuration failed.";
        return EXIT_FAILURE;
      }
    }

    // Export the used Image_describer and region type for:
    // - dynamic future regions computation and/or loading
    {
      std::ofstream stream(sImage_describer.c_str());
      if (!stream)
        return EXIT_FAILURE;

      cereal::JSONOutputArchive archive(stream);
      archive(cereal::make_nvp("image_describer", image_describer));
      auto regionsType = image_describer->Allocate();
      archive(cereal::make_nvp("regions_type", regionsType));
    }
  }

  // Feature extraction routines
  // For each View of the SfM_Data container:
  // - if regions file exists continue,
  // - if no file, compute features
  {
    system::Timer timer;
    Image<unsigned char> imageGray;

    system::LoggerProgress my_progress_bar(sfm_data.GetViews().size(), "- EXTRACT FEATURES -" );

    // Use a boolean to track if we must stop feature extraction
    std::atomic<bool> preemptive_exit(false);
  #ifdef OPENMVG_USE_OPENMP
      const unsigned int nb_max_thread = omp_get_max_threads();

      if (iNumThreads > 0) {
          omp_set_num_threads(iNumThreads);
      } else {
          omp_set_num_threads(nb_max_thread);
      }

      #pragma omp parallel for schedule(dynamic) if (iNumThreads > 0) private(imageGray)
  #endif
    for (int i = 0; i < static_cast<int>(sfm_data.views.size()); ++i)
    {
      Views::const_iterator iterViews = sfm_data.views.begin();
      std::advance(iterViews, i);
      const View * view = iterViews->second.get();
      const std::string
        sView_filename = stlplus::create_filespec(sfm_data.s_root_path, view->s_Img_path),
        sFeat = stlplus::create_filespec(sOutDir, stlplus::basename_part(sView_filename), "feat"),
        sDesc = stlplus::create_filespec(sOutDir, stlplus::basename_part(sView_filename), "desc");

      // If features or descriptors file are missing, compute them
      if (!preemptive_exit && (bForce || !stlplus::file_exists(sFeat) || !stlplus::file_exists(sDesc)))
      {
        if (!ReadImage(sView_filename.c_str(), &imageGray))
          continue;

        //
        // Look if there is an occlusion feature mask
        //
        Image<unsigned char> * mask = nullptr; // The mask is null by default

        const std::string
          mask_filename_local =
            stlplus::create_filespec(sfm_data.s_root_path,
              stlplus::basename_part(sView_filename) + "_mask", "png"),
          mask_filename_global =
            stlplus::create_filespec(sfm_data.s_root_path, "mask", "png");

        Image<unsigned char> imageMask;
        // Try to read the local mask
        if (stlplus::file_exists(mask_filename_local))
        {
          if (!ReadImage(mask_filename_local.c_str(), &imageMask))
          {
            OPENMVG_LOG_ERROR
              << "Invalid mask: " << mask_filename_local << ';'
              << "Stopping feature extraction.";
            preemptive_exit = true;
            continue;
          }
          // Use the local mask only if it fits the current image size
          if (imageMask.Width() == imageGray.Width() && imageMask.Height() == imageGray.Height())
            mask = &imageMask;
        }
        else
        {
          // Try to read the global mask
          if (stlplus::file_exists(mask_filename_global))
          {
            if (!ReadImage(mask_filename_global.c_str(), &imageMask))
            {
              OPENMVG_LOG_ERROR
                << "Invalid mask: " << mask_filename_global << ';'
                << "Stopping feature extraction.";
              preemptive_exit = true;
              continue;
            }
            // Use the global mask only if it fits the current image size
            if (imageMask.Width() == imageGray.Width() && imageMask.Height() == imageGray.Height())
              mask = &imageMask;
          }
        }

        // Compute features and descriptors and export them to files
        auto regions = image_describer->Describe(imageGray, mask);
        if (regions && !image_describer->Save(regions.get(), sFeat, sDesc)) {
          OPENMVG_LOG_ERROR
            << "Cannot save regions for image: " << sView_filename << ';'
            << "Stopping feature extraction.";
          preemptive_exit = true;
          continue;
        }
      }
      ++my_progress_bar;
    }
    OPENMVG_LOG_INFO << "Task done in (s): " << timer.elapsed();
  }
  return EXIT_SUCCESS;
}

// This executable computes pairs of images to be matched
int imagePairer(std::string sSfMDataFilename, std::string sOutputPairsFilename,
                std::string sPairMode, int iContiguousCount)
{
  // Mandatory elements:
  // cmd.add( make_option( 'i', sSfMDataFilename, "input_file" ) );
  // cmd.add( make_option( 'o', sOutputPairsFilename, "output_file" ) );
  // Optional elements:
  // cmd.add( make_option( 'm', sPairMode, "pair_mode" ) );
  // cmd.add( make_option( 'c', iContiguousCount, "contiguous_count" ) );

  // 0. Parse parameters
  std::cout << " You called:\n"
            // << argv[ 0 ] << "\n"
            << "--input_file       : " << sSfMDataFilename << "\n"
            << "--output_file      : " << sOutputPairsFilename << "\n"
            << "Optional parameters\n"
            << "--pair_mode        : " << sPairMode << "\n"
            << "--contiguous_count : " << iContiguousCount << "\n"
            << std::endl;

  EPairMode pairMode;
  if ( sPairMode == "EXHAUSTIVE" )
  {
    pairMode = PAIR_EXHAUSTIVE;
  }
  else if ( sPairMode == "CONTIGUOUS" )
  {
    if ( iContiguousCount == -1 )
    {
      // usage( argv[ 0 ] );
      std::cerr << "[Error] Contiguous pair mode selected but contiguous_count not set." << std::endl;
      exit( EXIT_FAILURE );
    }

    pairMode = PAIR_CONTIGUOUS;
  }

  // 1. Load SfM data scene
  std::cout << "Loading scene.";
  SfM_Data sfm_data;
  if ( !Load( sfm_data, sSfMDataFilename, ESfM_Data( VIEWS | INTRINSICS ) ) )
  {
    std::cerr << std::endl
              << "The input SfM_Data file \"" << sSfMDataFilename << "\" cannot be read." << std::endl;
    exit( EXIT_FAILURE );
  }
  const size_t NImage = sfm_data.GetViews().size();

  // 2. Compute pairs
  std::cout << "Computing pairs." << std::endl;
  Pair_Set pairs;
  switch ( pairMode )
  {
    case PAIR_EXHAUSTIVE:
    {
      pairs = exhaustivePairs( NImage );
      break;
    }
    case PAIR_CONTIGUOUS:
    {
      pairs = contiguousWithOverlap( NImage, iContiguousCount );
      break;
    }
    default:
    {
      std::cerr << "Unknown pair mode" << std::endl;
      exit( EXIT_FAILURE );
    }
  }

  // 3. Save pairs
  std::cout << "Saving pairs." << std::endl;
  if ( !savePairs( sOutputPairsFilename, pairs ) )
  {
    std::cerr << "Failed to save pairs to file: \"" << sOutputPairsFilename << "\"" << std::endl;
    exit( EXIT_FAILURE );
  }

  return EXIT_SUCCESS;
}

/// Compute corresponding features between a series of views:
/// - Load view images description (regions: features & descriptors)
/// - Compute putative local feature matches (descriptors matching)
int imageMatcher(std::string  sSfM_Data_Filename, std::string sOutputMatchesFilename, std::string  sPredefinedPairList,
                 float fDistRatio, std::string sNearestMatchingMethod, 
                 bool bForce, unsigned int ui_max_cache_size)
{
  // Pre-emptive matching parameters
  unsigned int ui_preemptive_feature_count = 200;
  double preemptive_matching_percentage_threshold = 0.08;
  bool pUsed = false;
  //required
  // cmd.add( make_option( 'i', sSfM_Data_Filename, "input_file" ) );
  // cmd.add( make_option( 'o', sOutputMatchesFilename, "output_file" ) );
  // cmd.add( make_option( 'p', sPredefinedPairList, "pair_list" ) );
  // // Options
  // cmd.add( make_option( 'r', fDistRatio, "ratio" ) );
  // cmd.add( make_option( 'n', sNearestMatchingMethod, "nearest_matching_method" ) );
  // cmd.add( make_option( 'f', bForce, "force" ) );
  // cmd.add( make_option( 'c', ui_max_cache_size, "cache_size" ) );
  // // Pre-emptive matching
  // cmd.add( make_option( 'P', ui_preemptive_feature_count, "preemptive_feature_count") );


  // try
  // {
  //   if ( argc == 1 )
  //     throw std::string( "Invalid command line parameter." );
  //   cmd.process( argc, argv );
  // }
  // catch ( const std::string& s )
  // {
  //   OPENMVG_LOG_INFO
  //     << "Usage: " << argv[ 0 ] << '\n'
  //     << "[-i|--input_file]   A SfM_Data file\n"
  //     << "[-o|--output_file]  Output file where computed matches are stored\n"
  //     << "[-p|--pair_list]    Pairs list file\n"
  //     << "\n[Optional]\n"
  //     << "[-f|--force] Force to recompute data]\n"
  //     << "[-r|--ratio] Distance ratio to discard non meaningful matches\n"
  //     << "   0.8: (default).\n"
  //     << "[-n|--nearest_matching_method]\n"
  //     << "  AUTO: auto choice from regions type,\n"
  //     << "  For Scalar based regions descriptor:\n"
  //     << "    BRUTEFORCEL2: L2 BruteForce matching,\n"
  //     << "    HNSWL2: L2 Approximate Matching with Hierarchical Navigable Small World graphs,\n"
  //     << "    HNSWL1: L1 Approximate Matching with Hierarchical Navigable Small World graphs\n"
  //     << "      tailored for quantized and histogram based descriptors (e.g uint8 RootSIFT)\n"
  //     << "    ANNL2: L2 Approximate Nearest Neighbor matching,\n"
  //     << "    CASCADEHASHINGL2: L2 Cascade Hashing matching.\n"
  //     << "    FASTCASCADEHASHINGL2: (default)\n"
  //     << "      L2 Cascade Hashing with precomputed hashed regions\n"
  //     << "     (faster than CASCADEHASHINGL2 but use more memory).\n"
  //     << "  For Binary based descriptor:\n"
  //     << "    BRUTEFORCEHAMMING: BruteForce Hamming matching,\n"
  //     << "    HNSWHAMMING: Hamming Approximate Matching with Hierarchical Navigable Small World graphs\n"
  //     << "[-c|--cache_size]\n"
  //     << "  Use a regions cache (only cache_size regions will be stored in memory)\n"
  //     << "  If not used, all regions will be load in memory."
  //     << "\n[Pre-emptive matching:]\n"
  //     << "[-P|--preemptive_feature_count] <NUMBER> Number of feature used for pre-emptive matching";

  //   OPENMVG_LOG_INFO << s;
  //   return EXIT_FAILURE;
  // }

  OPENMVG_LOG_INFO << " You called : "
            << "\n"
            // << argv[ 0 ] << "\n"
            << "--input_file " << sSfM_Data_Filename << "\n"
            << "--output_file " << sOutputMatchesFilename << "\n"
            << "--pair_list " << sPredefinedPairList << "\n"
            << "Optional parameters:"
            << "\n"
            << "--force " << bForce << "\n"
            << "--ratio " << fDistRatio << "\n"
            << "--nearest_matching_method " << sNearestMatchingMethod << "\n"
            << "--cache_size " << ((ui_max_cache_size == 0) ? "unlimited" : std::to_string(ui_max_cache_size)) << "\n"
            << "--preemptive_feature_used/count " << pUsed << " / " << ui_preemptive_feature_count;
  if (pUsed)
  {
    OPENMVG_LOG_INFO << "--preemptive_feature_count " << ui_preemptive_feature_count;
  }

  if ( sOutputMatchesFilename.empty() )
  {
    OPENMVG_LOG_ERROR << "No output file set.";
    return EXIT_FAILURE;
  }

  // -----------------------------
  // . Load SfM_Data Views & intrinsics data
  // . Compute putative descriptor matches
  // + Export some statistics
  // -----------------------------

  //---------------------------------------
  // Read SfM Scene (image view & intrinsics data)
  //---------------------------------------
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS))) {
    OPENMVG_LOG_ERROR << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read.";
    return EXIT_FAILURE;
  }
  const std::string sMatchesDirectory = stlplus::folder_part( sOutputMatchesFilename );

  //---------------------------------------
  // Load SfM Scene regions
  //---------------------------------------
  // Init the regions_type from the image describer file (used for image regions extraction)
  using namespace openMVG::features;
  const std::string sImage_describer = stlplus::create_filespec(sMatchesDirectory, "image_describer", "json");
  std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);
  if (!regions_type)
  {
    OPENMVG_LOG_ERROR << "Invalid: " << sImage_describer << " regions type file.";
    return EXIT_FAILURE;
  }

  //---------------------------------------
  // a. Compute putative descriptor matches
  //    - Descriptor matching (according user method choice)
  //    - Keep correspondences only if NearestNeighbor ratio is ok
  //---------------------------------------

  // Load the corresponding view regions
  std::shared_ptr<Regions_Provider> regions_provider;
  if (ui_max_cache_size == 0)
  {
    // Default regions provider (load & store all regions in memory)
    regions_provider = std::make_shared<Regions_Provider>();
  }
  else
  {
    // Cached regions provider (load & store regions on demand)
    regions_provider = std::make_shared<Regions_Provider_Cache>(ui_max_cache_size);
  }
  // If we use pre-emptive matching, we load less regions:
  if (ui_preemptive_feature_count > 0 && pUsed)
  {
    regions_provider = std::make_shared<Preemptive_Regions_Provider>(ui_preemptive_feature_count);
  }

  // Show the progress on the command line:
  system::LoggerProgress progress;

  if (!regions_provider->load(sfm_data, sMatchesDirectory, regions_type, &progress)) {
    OPENMVG_LOG_ERROR << "Cannot load view regions from: " << sMatchesDirectory << ".";
    return EXIT_FAILURE;
  }

  PairWiseMatches map_PutativeMatches;

  // Build some alias from SfM_Data Views data:
  // - List views as a vector of filenames & image sizes
  std::vector<std::string>               vec_fileNames;
  std::vector<std::pair<size_t, size_t>> vec_imagesSize;
  {
    vec_fileNames.reserve(sfm_data.GetViews().size());
    vec_imagesSize.reserve(sfm_data.GetViews().size());
    for (const auto view_it : sfm_data.GetViews())
    {
      const View * v = view_it.second.get();
      vec_fileNames.emplace_back(stlplus::create_filespec(sfm_data.s_root_path,
          v->s_Img_path));
      vec_imagesSize.emplace_back(v->ui_width, v->ui_height);
    }
  }

  OPENMVG_LOG_INFO << " - PUTATIVE MATCHES - ";
  // If the matches already exists, reload them
  if ( !bForce && ( stlplus::file_exists( sOutputMatchesFilename ) ) )
  {
    if ( !( Load( map_PutativeMatches, sOutputMatchesFilename ) ) )
    {
      OPENMVG_LOG_ERROR << "Cannot load input matches file";
      return EXIT_FAILURE;
    }
    OPENMVG_LOG_INFO
      << "\t PREVIOUS RESULTS LOADED;"
      << " #pair: " << map_PutativeMatches.size();
  }
  else // Compute the putative matches
  {
    // Allocate the right Matcher according the Matching requested method
    std::unique_ptr<Matcher> collectionMatcher;
    if ( sNearestMatchingMethod == "AUTO" )
    {
      if ( regions_type->IsScalar() )
      {
        OPENMVG_LOG_INFO << "Using FAST_CASCADE_HASHING_L2 matcher";
        collectionMatcher.reset(new Cascade_Hashing_Matcher_Regions(fDistRatio));
      }
      else
      if (regions_type->IsBinary())
      {
        OPENMVG_LOG_INFO << "Using HNSWHAMMING matcher";
        collectionMatcher.reset(new Matcher_Regions(fDistRatio, HNSW_HAMMING));
      }
    }
    else
    if (sNearestMatchingMethod == "BRUTEFORCEL2")
    {
      OPENMVG_LOG_INFO << "Using BRUTE_FORCE_L2 matcher";
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, BRUTE_FORCE_L2));
    }
    else
    if (sNearestMatchingMethod == "BRUTEFORCEHAMMING")
    {
      OPENMVG_LOG_INFO << "Using BRUTE_FORCE_HAMMING matcher";
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, BRUTE_FORCE_HAMMING));
    }
    else
    if (sNearestMatchingMethod == "HNSWL2")
    {
      OPENMVG_LOG_INFO << "Using HNSWL2 matcher";
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, HNSW_L2));
    }
    if (sNearestMatchingMethod == "HNSWL1")
    {
      OPENMVG_LOG_INFO << "Using HNSWL1 matcher";
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, HNSW_L1));
    }
    else
    if (sNearestMatchingMethod == "HNSWHAMMING")
    {
      OPENMVG_LOG_INFO << "Using HNSWHAMMING matcher";
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, HNSW_HAMMING));
    }
    else
    if (sNearestMatchingMethod == "ANNL2")
    {
      OPENMVG_LOG_INFO << "Using ANN_L2 matcher";
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, ANN_L2));
    }
    else
    if (sNearestMatchingMethod == "CASCADEHASHINGL2")
    {
      OPENMVG_LOG_INFO << "Using CASCADE_HASHING_L2 matcher";
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, CASCADE_HASHING_L2));
    }
    else
    if (sNearestMatchingMethod == "FASTCASCADEHASHINGL2")
    {
      OPENMVG_LOG_INFO << "Using FAST_CASCADE_HASHING_L2 matcher";
      collectionMatcher.reset(new Cascade_Hashing_Matcher_Regions(fDistRatio));
    }
    if (!collectionMatcher)
    {
      OPENMVG_LOG_ERROR << "Invalid Nearest Neighbor method: " << sNearestMatchingMethod;
      return EXIT_FAILURE;
    }
    // Perform the matching
    system::Timer timer;
    {
      // From matching mode compute the pair list that have to be matched:
      Pair_Set pairs;
      if ( sPredefinedPairList.empty() )
      {
        OPENMVG_LOG_INFO << "No input pair file set. Use exhaustive match by default.";
        const size_t NImage = sfm_data.GetViews().size();
        pairs = exhaustivePairs( NImage );
      }
      else
      if ( !loadPairs( sfm_data.GetViews().size(), sPredefinedPairList, pairs ) )
      {
        OPENMVG_LOG_ERROR << "Failed to load pairs from file: \"" << sPredefinedPairList << "\"";
        return EXIT_FAILURE;
      }
      OPENMVG_LOG_INFO << "Running matching on #pairs: " << pairs.size();
      // Photometric matching of putative pairs
      OPENMVG_LOG_INFO << "this line?";
      collectionMatcher->Match( regions_provider, pairs, map_PutativeMatches, &progress );
      OPENMVG_LOG_INFO << pUsed;
      if (pUsed) // Preemptive filter
      {
        // Keep putative matches only if there is more than X matches
        PairWiseMatches map_filtered_matches;
        for (const auto & pairwisematches_it : map_PutativeMatches)
        {
          const size_t putative_match_count = pairwisematches_it.second.size();
          const int match_count_threshold =
            preemptive_matching_percentage_threshold * ui_preemptive_feature_count;
          // TODO: Add an option to keeping X Best pairs
          if (putative_match_count >= match_count_threshold)  {
            // the pair will be kept
            map_filtered_matches.insert(pairwisematches_it);
          }
        }
        map_PutativeMatches.clear();
        std::swap(map_filtered_matches, map_PutativeMatches);
      }
      OPENMVG_LOG_INFO << "put?";
      //---------------------------------------
      //-- Export putative matches & pairs
      //---------------------------------------
      if ( !Save( map_PutativeMatches, std::string( sOutputMatchesFilename ) ) )
      {
        OPENMVG_LOG_ERROR
          << "Cannot save computed matches in: "
          << sOutputMatchesFilename;
        return EXIT_FAILURE;
      }

      OPENMVG_LOG_INFO << "save?";
      // Save pairs
      const std::string sOutputPairFilename =
        stlplus::create_filespec( sMatchesDirectory, "preemptive_pairs", "txt" );
      if (!savePairs(
        sOutputPairFilename,
        getPairs(map_PutativeMatches)))
      {
        OPENMVG_LOG_ERROR
          << "Cannot save computed matches pairs in: "
          << sOutputPairFilename;
        return EXIT_FAILURE;
      }
      OPENMVG_LOG_INFO << "end region";
    }
    OPENMVG_LOG_INFO << "Task (Regions Matching) done in (s): " << timer.elapsed();
  }

  OPENMVG_LOG_INFO << "#Putative pairs: " << map_PutativeMatches.size();

  // -- export Putative View Graph statistics
  graph::getGraphStatistics(sfm_data.GetViews().size(), getPairs(map_PutativeMatches));

  //-- export putative matches Adjacency matrix
  PairWiseMatchingToAdjacencyMatrixSVG( vec_fileNames.size(),
                                        map_PutativeMatches,
                                        stlplus::create_filespec( sMatchesDirectory, "PutativeAdjacencyMatrix", "svg" ) );
  //-- export view pair graph once putative graph matches has been computed
  {
    std::set<IndexT> set_ViewIds;
    std::transform( sfm_data.GetViews().begin(), sfm_data.GetViews().end(), std::inserter( set_ViewIds, set_ViewIds.begin() ), stl::RetrieveKey() );
    graph::indexedGraph putativeGraph( set_ViewIds, getPairs( map_PutativeMatches ) );
    graph::exportToGraphvizData(
        stlplus::create_filespec( sMatchesDirectory, "putative_matches" ),
        putativeGraph );
  }

  return EXIT_SUCCESS;
}

enum EGeometricModel
{
  FUNDAMENTAL_MATRIX       = 0,
  ESSENTIAL_MATRIX         = 1,
  HOMOGRAPHY_MATRIX        = 2,
  ESSENTIAL_MATRIX_ANGULAR = 3,
  ESSENTIAL_MATRIX_ORTHO   = 4,
  ESSENTIAL_MATRIX_UPRIGHT = 5
};

/// Compute corresponding features between a series of views:
/// - Load view images description (regions: features & descriptors)
/// - Compute putative local feature matches (descriptors matching)
/// - Compute geometric coherent feature matches (robust model estimation from putative matches)
/// - Export computed data
int imageFilterGeo(std::string sSfM_Data_Filename, std::string sPutativeMatchesFilename, 
                    std::string sFilteredMatchesFilename)
{
  // // The scene
  // std::string sSfM_Data_Filename;
  // // The input matches
  // std::string sPutativeMatchesFilename;
  // // The output matches
  // std::string sFilteredMatchesFilename;
  // The input pairs
  std::string sInputPairsFilename;
  // The output pairs
  std::string sOutputPairsFilename;

  std::string  sGeometricModel   = "f";
  bool         bForce            = false;
  bool         bGuided_matching  = false;
  int          imax_iteration    = 2048;
  unsigned int ui_max_cache_size = 0;

  //required
  // cmd.add( make_option( 'i', sSfM_Data_Filename, "input_file" ) );
  // cmd.add( make_option( 'o', sFilteredMatchesFilename, "output_file" ) );
  // cmd.add( make_option( 'm', sPutativeMatchesFilename, "matches" ) );
  // // Options
  // cmd.add( make_option( 'p', sInputPairsFilename, "input_pairs" ) );
  // cmd.add( make_option( 's', sOutputPairsFilename, "output_pairs" ) );
  // cmd.add( make_option( 'g', sGeometricModel, "geometric_model" ) );
  // cmd.add( make_option( 'f', bForce, "force" ) );
  // cmd.add( make_option( 'r', bGuided_matching, "guided_matching" ) );
  // cmd.add( make_option( 'I', imax_iteration, "max_iteration" ) );
  // cmd.add( make_option( 'c', ui_max_cache_size, "cache_size" ) );

  // try
  // {
  //   if ( argc == 1 )
  //     throw std::string( "Invalid command line parameter." );
  //   cmd.process( argc, argv );
  // }
  // catch ( const std::string& s )
  // {
  //   OPENMVG_LOG_INFO << "Usage: " << argv[0] << '\n'
  //                    << "[-i|--input_file]       A SfM_Data file\n"
  //                    << "[-m|--matches]          (Input) matches filename\n"
  //                    << "[-o|--output_file]      (Output) filtered matches filename\n"
  //                    << "\n[Optional]\n"
  //                    << "[-p|--input_pairs]      (Input) pairs filename\n"
  //                    << "[-s|--output_pairs]     (Output) filtered pairs filename\n"
  //                    << "[-f|--force]            Force to recompute data\n"
  //                    << "[-g|--geometric_model]\n"
  //                    << "  (pairwise correspondences filtering thanks to robust model estimation):\n"
  //                    << "   f: (default) fundamental matrix,\n"
  //                    << "   e: essential matrix,\n"
  //                    << "   h: homography matrix.\n"
  //                    << "   a: essential matrix with an angular parametrization,\n"
  //                    << "   u: upright essential matrix with an angular parametrization,\n"
  //                    << "   o: orthographic essential matrix.\n"
  //                    << "[-r|--guided_matching]  Use the found model to improve the pairwise correspondences.\n"
  //                    << "[-c|--cache_size]\n"
  //                    << "  Use a regions cache (only cache_size regions will be stored in memory)\n"
  //                    << "  If not used, all regions will be load in memory.";

  //   OPENMVG_LOG_INFO << s;
  //   return EXIT_FAILURE;
  // }

  OPENMVG_LOG_INFO << " You called : "
                   << "\n"
                  //  << argv[0] << "\n"
                   << "--input_file:        " << sSfM_Data_Filename << "\n"
                   << "--matches:           " << sPutativeMatchesFilename << "\n"
                   << "--output_file:       " << sFilteredMatchesFilename << "\n"
                   << "Optional parameters: "
                   << "\n"
                   << "--input_pairs        " << sInputPairsFilename << "\n"
                   << "--output_pairs       " << sOutputPairsFilename << "\n"
                   << "--force              " << (bForce ? "true" : "false") << "\n"
                   << "--geometric_model    " << sGeometricModel << "\n"
                   << "--guided_matching    " << bGuided_matching << "\n"
                   << "--cache_size         " << ((ui_max_cache_size == 0) ? "unlimited" : std::to_string(ui_max_cache_size));

  if ( sFilteredMatchesFilename.empty() )
  {
    OPENMVG_LOG_ERROR << "It is an invalid output file";
    return EXIT_FAILURE;
  }
  if ( sSfM_Data_Filename.empty() )
  {
    OPENMVG_LOG_ERROR << "It is an invalid SfM file";
    return EXIT_FAILURE;
  }
  if ( sPutativeMatchesFilename.empty() )
  {
    OPENMVG_LOG_ERROR << "It is an invalid putative matche file";
    return EXIT_FAILURE;
  }

  const std::string sMatchesDirectory = stlplus::folder_part( sPutativeMatchesFilename );

  EGeometricModel eGeometricModelToCompute = FUNDAMENTAL_MATRIX;
  switch ( std::tolower(sGeometricModel[ 0 ], std::locale()) )
  {
    case 'f':
      eGeometricModelToCompute = FUNDAMENTAL_MATRIX;
      break;
    case 'e':
      eGeometricModelToCompute = ESSENTIAL_MATRIX;
      break;
    case 'h':
      eGeometricModelToCompute = HOMOGRAPHY_MATRIX;
      break;
    case 'a':
      eGeometricModelToCompute = ESSENTIAL_MATRIX_ANGULAR;
      break;
    case 'u':
      eGeometricModelToCompute = ESSENTIAL_MATRIX_UPRIGHT;
      break;
    case 'o':
      eGeometricModelToCompute = ESSENTIAL_MATRIX_ORTHO;
      break;
    default:
      OPENMVG_LOG_ERROR << "Unknown geometric model";
      return EXIT_FAILURE;
  }

  // -----------------------------
  // - Load SfM_Data Views & intrinsics data
  // a. Load putative descriptor matches
  // [a.1] Filter matches with input pairs
  // b. Geometric filtering of putative matches
  // + Export some statistics
  // -----------------------------

  //---------------------------------------
  // Read SfM Scene (image view & intrinsics data)
  //---------------------------------------
  SfM_Data sfm_data;
  if ( !Load( sfm_data, sSfM_Data_Filename, ESfM_Data( VIEWS | INTRINSICS ) ) )
  {
    OPENMVG_LOG_ERROR << "The input SfM_Data file \"" << sSfM_Data_Filename << "\" cannot be read.";
    return EXIT_FAILURE;
  }

  //---------------------------------------
  // Load SfM Scene regions
  //---------------------------------------
  // Init the regions_type from the image describer file (used for image regions extraction)
  using namespace openMVG::features;
  // Consider that the image_describer.json is inside the matches directory (which is bellow the sfm_data.bin)
  const std::string        sImage_describer = stlplus::create_filespec( sMatchesDirectory, "image_describer.json" );
  std::unique_ptr<Regions> regions_type     = Init_region_type_from_file( sImage_describer );
  if ( !regions_type )
  {
    OPENMVG_LOG_ERROR << "Invalid: " << sImage_describer << " regions type file.";
    return EXIT_FAILURE;
  }

  //---------------------------------------
  // a. Compute putative descriptor matches
  //    - Descriptor matching (according user method choice)
  //    - Keep correspondences only if NearestNeighbor ratio is ok
  //---------------------------------------

  // Load the corresponding view regions
  std::shared_ptr<Regions_Provider> regions_provider;
  if ( ui_max_cache_size == 0 )
  {
    // Default regions provider (load & store all regions in memory)
    regions_provider = std::make_shared<Regions_Provider>();
  }
  else
  {
    // Cached regions provider (load & store regions on demand)
    regions_provider = std::make_shared<Regions_Provider_Cache>( ui_max_cache_size );
  }

  // Show the progress on the command line:
  system::LoggerProgress progress;

  if ( !regions_provider->load( sfm_data, sMatchesDirectory, regions_type, &progress ) )
  {
    OPENMVG_LOG_ERROR << "Invalid regions.";
    return EXIT_FAILURE;
  }

  PairWiseMatches map_PutativeMatches;
  //---------------------------------------
  // A. Load initial matches
  //---------------------------------------
  if ( !Load( map_PutativeMatches, sPutativeMatchesFilename ) )
  {
    OPENMVG_LOG_ERROR << "Failed to load the initial matches file.";
    return EXIT_FAILURE;
  }

  if ( !sInputPairsFilename.empty() )
  {
    // Load input pairs
    OPENMVG_LOG_INFO << "Loading input pairs ...";
    Pair_Set input_pairs;
    loadPairs( sfm_data.GetViews().size(), sInputPairsFilename, input_pairs );

    // Filter matches with the given pairs
    OPENMVG_LOG_INFO << "Filtering matches with the given pairs.";
    map_PutativeMatches = getPairs( map_PutativeMatches, input_pairs );
  }

  //---------------------------------------
  // b. Geometric filtering of putative matches
  //    - AContrario Estimation of the desired geometric model
  //    - Use an upper bound for the a contrario estimated threshold
  //---------------------------------------

  std::unique_ptr<ImageCollectionGeometricFilter> filter_ptr(
      new ImageCollectionGeometricFilter( &sfm_data, regions_provider ) );

  if ( filter_ptr )
  {
    system::Timer timer;
    const double  d_distance_ratio = 0.6;

    PairWiseMatches map_GeometricMatches;
    switch ( eGeometricModelToCompute )
    {
      case HOMOGRAPHY_MATRIX:
      {
        const bool bGeometric_only_guided_matching = true;
        filter_ptr->Robust_model_estimation(
            GeometricFilter_HMatrix_AC( 4.0, imax_iteration ),
            map_PutativeMatches,
            bGuided_matching,
            bGeometric_only_guided_matching ? -1.0 : d_distance_ratio,
            &progress );
        map_GeometricMatches = filter_ptr->Get_geometric_matches();
      }
      break;
      case FUNDAMENTAL_MATRIX:
      {
        filter_ptr->Robust_model_estimation(
            GeometricFilter_FMatrix_AC( 4.0, imax_iteration ),
            map_PutativeMatches,
            bGuided_matching,
            d_distance_ratio,
            &progress );
        map_GeometricMatches = filter_ptr->Get_geometric_matches();
      }
      break;
      case ESSENTIAL_MATRIX:
      {
        filter_ptr->Robust_model_estimation(
            GeometricFilter_EMatrix_AC( 4.0, imax_iteration ),
            map_PutativeMatches,
            bGuided_matching,
            d_distance_ratio,
            &progress );
        map_GeometricMatches = filter_ptr->Get_geometric_matches();

        //-- Perform an additional check to remove pairs with poor overlap
        std::vector<PairWiseMatches::key_type> vec_toRemove;
        for ( const auto& pairwisematches_it : map_GeometricMatches )
        {
          const size_t putativePhotometricCount = map_PutativeMatches.find( pairwisematches_it.first )->second.size();
          const size_t putativeGeometricCount   = pairwisematches_it.second.size();
          const float  ratio                    = putativeGeometricCount / static_cast<float>( putativePhotometricCount );
          if ( putativeGeometricCount < 50 || ratio < .3f )
          {
            // the pair will be removed
            vec_toRemove.push_back( pairwisematches_it.first );
          }
        }
        //-- remove discarded pairs
        for ( const auto& pair_to_remove_it : vec_toRemove )
        {
          map_GeometricMatches.erase( pair_to_remove_it );
        }
      }
      break;
      case ESSENTIAL_MATRIX_ANGULAR:
      {
        filter_ptr->Robust_model_estimation(
          GeometricFilter_ESphericalMatrix_AC_Angular<false>(4.0, imax_iteration),
          map_PutativeMatches, bGuided_matching, d_distance_ratio, &progress);
        map_GeometricMatches = filter_ptr->Get_geometric_matches();
      }
      break;
      case ESSENTIAL_MATRIX_UPRIGHT:
      {
        filter_ptr->Robust_model_estimation(
          GeometricFilter_ESphericalMatrix_AC_Angular<true>(4.0, imax_iteration),
          map_PutativeMatches, bGuided_matching, d_distance_ratio, &progress);
        map_GeometricMatches = filter_ptr->Get_geometric_matches();
      }
      break;
      case ESSENTIAL_MATRIX_ORTHO:
      {
        filter_ptr->Robust_model_estimation(
            GeometricFilter_EOMatrix_RA( 2.0, imax_iteration ),
            map_PutativeMatches,
            bGuided_matching,
            d_distance_ratio,
            &progress );
        map_GeometricMatches = filter_ptr->Get_geometric_matches();
      }
      break;
    }

    //---------------------------------------
    //-- Export geometric filtered matches
    //---------------------------------------
    if ( !Save( map_GeometricMatches, sFilteredMatchesFilename ) )
    {
      OPENMVG_LOG_ERROR << "Cannot save filtered matches in: " << sFilteredMatchesFilename;
      return EXIT_FAILURE;
    }

    // -- export Geometric View Graph statistics
    graph::getGraphStatistics(sfm_data.GetViews().size(), getPairs(map_GeometricMatches));

    OPENMVG_LOG_INFO << "Task done in (s): " << timer.elapsed();

    //-- export Adjacency matrix
    OPENMVG_LOG_INFO <<  "\n Export Adjacency Matrix of the pairwise's geometric matches";

    PairWiseMatchingToAdjacencyMatrixSVG( sfm_data.GetViews().size(),
                                          map_GeometricMatches,
                                          stlplus::create_filespec( sMatchesDirectory, "GeometricAdjacencyMatrix", "svg" ) );

    const Pair_Set outputPairs = getPairs( map_GeometricMatches );

    //-- export view pair graph once geometric filter have been done
    {
      std::set<IndexT> set_ViewIds;
      std::transform( sfm_data.GetViews().begin(), sfm_data.GetViews().end(), std::inserter( set_ViewIds, set_ViewIds.begin() ), stl::RetrieveKey() );
      graph::indexedGraph putativeGraph( set_ViewIds, outputPairs );
      graph::exportToGraphvizData(
          stlplus::create_filespec( sMatchesDirectory, "geometric_matches" ),
          putativeGraph );
    }

    // Write pairs
    if ( !sOutputPairsFilename.empty() )
    {
      OPENMVG_LOG_INFO << "Saving pairs to: " << sOutputPairsFilename;
      if ( !savePairs( sOutputPairsFilename, outputPairs ) )
      {
        OPENMVG_LOG_ERROR << "Failed to write pairs file";
        return EXIT_FAILURE;
      }
    }
  }
  return EXIT_SUCCESS;
}

enum class ESfMSceneInitializer
{
  INITIALIZE_EXISTING_POSES,
  INITIALIZE_MAX_PAIR,
  INITIALIZE_AUTO_PAIR,
  INITIALIZE_STELLAR
};

enum class ESfMEngine
{
  INCREMENTAL,
  INCREMENTALV2,
  GLOBAL
};

bool StringToEnum
(
  const std::string & str,
  ESfMEngine & sfm_engine
)
{
  const std::map<std::string, ESfMEngine> string_to_enum_mapping =
  {
    {"INCREMENTAL", ESfMEngine::INCREMENTAL},
    {"INCREMENTALV2", ESfMEngine::INCREMENTALV2},
    {"GLOBAL", ESfMEngine::GLOBAL},
  };
  const auto it  = string_to_enum_mapping.find(str);
  if (it == string_to_enum_mapping.end())
    return false;
  sfm_engine = it->second;
  return true;
}

bool StringToEnum
(
  const std::string & str,
  ESfMSceneInitializer & scene_initializer
)
{
  const std::map<std::string, ESfMSceneInitializer> string_to_enum_mapping =
  {
    {"EXISTING_POSE", ESfMSceneInitializer::INITIALIZE_EXISTING_POSES},
    {"MAX_PAIR", ESfMSceneInitializer::INITIALIZE_MAX_PAIR},
    {"AUTO_PAIR", ESfMSceneInitializer::INITIALIZE_AUTO_PAIR},
    {"STELLAR", ESfMSceneInitializer::INITIALIZE_STELLAR},
  };
  const auto it  = string_to_enum_mapping.find(str);
  if (it == string_to_enum_mapping.end())
    return false;
  scene_initializer = it->second;
  return true;
}

/// From 2 given image filenames, find the two corresponding index in the View list
bool computeIndexFromImageNames(
  const SfM_Data & sfm_data,
  const std::pair<std::string,std::string>& initialPairName,
  Pair& initialPairIndex)
{
  if (initialPairName.first == initialPairName.second)
  {
    OPENMVG_LOG_ERROR << "Invalid image names. You cannot use the same image to initialize a pair.";
    return false;
  }

  initialPairIndex = {UndefinedIndexT, UndefinedIndexT};

  /// List views filenames and find the one that correspond to the user ones:
  for (Views::const_iterator it = sfm_data.GetViews().begin();
     it != sfm_data.GetViews().end(); ++it)
  {
    const View * v = it->second.get();
    const std::string filename = stlplus::filename_part(v->s_Img_path);
    if (filename == initialPairName.first)
    {
      initialPairIndex.first = v->id_view;
    }
    else
    {
      if (filename == initialPairName.second)
      {
        initialPairIndex.second = v->id_view;
      }
    }
  }
  return (initialPairIndex.first != UndefinedIndexT &&
      initialPairIndex.second != UndefinedIndexT);
}

int imageSfm(std::string filename_sfm_data, std::string directory_match, std::string directory_output,
             std::string engine_name)
{
  // OPENMVG_LOG_INFO
  //     << "\n-----------------------------------------------------------"
  //     << "\n Structure from Motion:"
  //     << "\n-----------------------------------------------------------";
  // CmdLine cmd;

  // Common options:
  // std::string
  //     filename_sfm_data,
  //     directory_match,
    std::string  filename_match;
      // directory_output,
      // engine_name = "INCREMENTAL";

  // Bundle adjustment options:
  std::string sIntrinsic_refinement_options = "ADJUST_ALL";
  std::string sExtrinsic_refinement_options = "ADJUST_ALL";
  bool b_use_motion_priors = false;

  // Incremental SfM options
  int triangulation_method = static_cast<int>(ETriangulationMethod::DEFAULT);
  int resection_method  = static_cast<int>(resection::SolverType::DEFAULT);
  int user_camera_model = PINHOLE_CAMERA_RADIAL3;

  // SfM v1
  std::pair<std::string,std::string> initial_pair_string("","");

  // SfM v2
  std::string sfm_initializer_method = "STELLAR";

  // Global SfM
  int rotation_averaging_method = int (ROTATION_AVERAGING_L2);
  int translation_averaging_method = int (TRANSLATION_AVERAGING_SOFTL1);


  // Common options
  // cmd.add( make_option('i', filename_sfm_data, "input_file") );
  // cmd.add( make_option('m', directory_match, "match_dir") );
  // cmd.add( make_option('M', filename_match, "match_file") );
  // cmd.add( make_option('o', directory_output, "output_dir") );
  // cmd.add( make_option('s', engine_name, "sfm_engine") );

  // // Bundle adjustment options
  // cmd.add( make_option('f', sIntrinsic_refinement_options, "refine_intrinsic_config") );
  // cmd.add( make_option('e', sExtrinsic_refinement_options, "refine_extrinsic_config") );
  // cmd.add( make_switch('P', "prior_usage") );

  // // Incremental SfM pipeline options
  // cmd.add( make_option('t', triangulation_method, "triangulation_method"));
  // cmd.add( make_option('r', resection_method, "resection_method"));
  // cmd.add( make_option('c', user_camera_model, "camera_model") );
  // // Incremental SfM2
  // cmd.add( make_option('S', sfm_initializer_method, "sfm_initializer") );
  // // Incremental SfM1
  // cmd.add( make_option('a', initial_pair_string.first, "initial_pair_a") );
  // cmd.add( make_option('b', initial_pair_string.second, "initial_pair_b") );
  // // Global SfM
  // cmd.add( make_option('R', rotation_averaging_method, "rotationAveraging") );
  // cmd.add( make_option('T', translation_averaging_method, "translationAveraging") );

  // try {
  //   if (argc == 1) throw std::string("Invalid parameter.");
  //   cmd.process(argc, argv);
  // } catch (const std::string& s) {

  //   OPENMVG_LOG_INFO << "Usage: " << argv[0] << '\n'
  //     << "[Required]\n"
  //     << "[-i|--input_file] path to a SfM_Data scene\n"
  //     << "[-m|--match_dir] path to the matches that corresponds to the provided SfM_Data scene\n"
  //     << "[-o|--output_dir] path where the output data will be stored\n"
  //     << "[-s|--sfm_engine] Type of SfM Engine to use for the reconstruction\n"
  //     << "\t INCREMENTAL   : add image sequentially to a 2 view seed\n"
  //     << "\t INCREMENTALV2 : add image sequentially to a 2 or N view seed (experimental)\n"
  //     << "\t GLOBAL    : initialize globally rotation and translations\n"
  //     << "\n\n"
  //     << "[Optional parameters]\n"
  //     << "\n\n"
  //     << "[Common]\n"
  //     << "[-M|--match_file] path to the match file to use (i.e matches.f.txt or matches.f.bin)\n"
  //     << "[-f|--refine_intrinsic_config] Intrinsic parameters refinement option\n"
  //     << "\t ADJUST_ALL -> refine all existing parameters (default) \n"
  //     << "\t NONE -> intrinsic parameters are held as constant\n"
  //     << "\t ADJUST_FOCAL_LENGTH -> refine only the focal length\n"
  //     << "\t ADJUST_PRINCIPAL_POINT -> refine only the principal point position\n"
  //     << "\t ADJUST_DISTORTION -> refine only the distortion coefficient(s) (if any)\n"
  //     << "\t -> NOTE: options can be combined thanks to '|'\n"
  //     << "\t ADJUST_FOCAL_LENGTH|ADJUST_PRINCIPAL_POINT\n"
  //     <<    "\t\t-> refine the focal length & the principal point position\n"
  //     << "\t ADJUST_FOCAL_LENGTH|ADJUST_DISTORTION\n"
  //     <<    "\t\t-> refine the focal length & the distortion coefficient(s) (if any)\n"
  //     << "\t ADJUST_PRINCIPAL_POINT|ADJUST_DISTORTION\n"
  //     <<    "\t\t-> refine the principal point position & the distortion coefficient(s) (if any)\n"
  //     << "[-e|--refine_extrinsic_config] Extrinsic parameters refinement option\n"
  //     << "\t ADJUST_ALL -> refine all existing parameters (default) \n"
  //     << "\t NONE -> extrinsic parameters are held as constant\n"
  //     << "[-P|--prior_usage] Enable usage of motion priors (i.e GPS positions) (default: false)\n"
  //     << "\n\n"
  //     << "[Engine specifics]\n"
  //     << "\n\n"
  //     << "[INCREMENTAL]\n"
  //     << "\t[-a|--initial_pair_a] filename of the first image (without path)\n"
  //     << "\t[-b|--initial_pair_b] filename of the second image (without path)\n"
  //     << "\t[-c|--camera_model] Camera model type for view with unknown intrinsic:\n"
  //     << "\t\t 1: Pinhole \n"
  //     << "\t\t 2: Pinhole radial 1\n"
  //     << "\t\t 3: Pinhole radial 3 (default)\n"
  //     << "\t\t 4: Pinhole radial 3 + tangential 2\n"
  //     << "\t\t 5: Pinhole fisheye\n"
  //     << "\t[--triangulation_method] triangulation method (default=" << triangulation_method << "):\n"
  //     << "\t\t" << static_cast<int>(ETriangulationMethod::DIRECT_LINEAR_TRANSFORM) << ": DIRECT_LINEAR_TRANSFORM\n"
  //     << "\t\t" << static_cast<int>(ETriangulationMethod::L1_ANGULAR) << ": L1_ANGULAR\n"
  //     << "\t\t" << static_cast<int>(ETriangulationMethod::LINFINITY_ANGULAR) << ": LINFINITY_ANGULAR\n"
  //     << "\t\t" << static_cast<int>(ETriangulationMethod::INVERSE_DEPTH_WEIGHTED_MIDPOINT) << ": INVERSE_DEPTH_WEIGHTED_MIDPOINT\n"
  //     << "\t[--resection_method] resection/pose estimation method (default=" << resection_method << "):\n"
  //     << "\t\t" << static_cast<int>(resection::SolverType::DLT_6POINTS) << ": DIRECT_LINEAR_TRANSFORM 6Points | does not use intrinsic data\n"
  //     << "\t\t" << static_cast<int>(resection::SolverType::P3P_KE_CVPR17) << ": P3P_KE_CVPR17\n"
  //     << "\t\t" << static_cast<int>(resection::SolverType::P3P_KNEIP_CVPR11) << ": P3P_KNEIP_CVPR11\n"
  //     << "\t\t" << static_cast<int>(resection::SolverType::P3P_NORDBERG_ECCV18) << ": P3P_NORDBERG_ECCV18\n"
  //     << "\t\t" << static_cast<int>(resection::SolverType::UP2P_KUKELOVA_ACCV10)  << ": UP2P_KUKELOVA_ACCV10 | 2Points | upright camera\n"
  //     << "\n\n"
  //     << "[INCREMENTALV2]\n"
  //     << "\t[-S|--sfm_initializer] Choose the SfM initializer method:\n"
  //     << "\t\t 'EXISTING_POSE'-> Initialize the reconstruction from the existing sfm_data camera poses\n"
  //     << "\t\t 'MAX_PAIR'-> Initialize the reconstruction from the pair that has the most of matches\n"
  //     << "\t\t 'AUTO_PAIR'-> Initialize the reconstruction with a pair selected automatically\n"
  //     << "\t\t 'STELLAR'-> Initialize the reconstruction with a 'stellar' reconstruction\n"
  //     << "\t[-c|--camera_model] Camera model type for view with unknown intrinsic:\n"
  //     << "\t\t 1: Pinhole \n"
  //     << "\t\t 2: Pinhole radial 1\n"
  //     << "\t\t 3: Pinhole radial 3 (default)\n"
  //     << "\t\t 4: Pinhole radial 3 + tangential 2\n"
  //     << "\t\t 5: Pinhole fisheye\n"
  //     << "\t[--triangulation_method] triangulation method (default=" << triangulation_method << "):\n"
  //     << "\t\t" << static_cast<int>(ETriangulationMethod::DIRECT_LINEAR_TRANSFORM) << ": DIRECT_LINEAR_TRANSFORM\n"
  //     << "\t\t" << static_cast<int>(ETriangulationMethod::L1_ANGULAR) << ": L1_ANGULAR\n"
  //     << "\t\t" << static_cast<int>(ETriangulationMethod::LINFINITY_ANGULAR) << ": LINFINITY_ANGULAR\n"
  //     << "\t\t" << static_cast<int>(ETriangulationMethod::INVERSE_DEPTH_WEIGHTED_MIDPOINT) << ": INVERSE_DEPTH_WEIGHTED_MIDPOINT\n"
  //     << "\t[--resection_method] resection/pose estimation method (default=" << resection_method << "):\n"
  //     << "\t\t" << static_cast<int>(resection::SolverType::DLT_6POINTS) << ": DIRECT_LINEAR_TRANSFORM 6Points | does not use intrinsic data\n"
  //     << "\t\t" << static_cast<int>(resection::SolverType::P3P_KE_CVPR17) << ": P3P_KE_CVPR17\n"
  //     << "\t\t" << static_cast<int>(resection::SolverType::P3P_KNEIP_CVPR11) << ": P3P_KNEIP_CVPR11\n"
  //     << "\t\t" << static_cast<int>(resection::SolverType::P3P_NORDBERG_ECCV18) << ": P3P_NORDBERG_ECCV18\n"
  //     << "\t\t" << static_cast<int>(resection::SolverType::UP2P_KUKELOVA_ACCV10)  << ": UP2P_KUKELOVA_ACCV10 | 2Points | upright camera\n"
  //     << "\n\n"
  //     << "[GLOBAL]\n"
  //     << "\t[-R|--rotationAveraging]\n"
  //     << "\t\t 1 -> L1 minimization\n"
  //     << "\t\t 2 -> L2 minimization (default)\n"
  //     << "\t[-T|--translationAveraging]:\n"
  //     << "\t\t 1 -> L1 minimization\n"
  //     << "\t\t 2 -> L2 minimization of sum of squared Chordal distances\n"
  //     << "\t\t 3 -> SoftL1 minimization (default)\n"
  //     << "\t\t 4 -> LiGT: Linear Global Translation constraints from rotation and matches\n";

  //   OPENMVG_LOG_ERROR << s;
  //   return EXIT_FAILURE;
  // }

  // b_use_motion_priors = cmd.used('P');
  b_use_motion_priors = false;

  // Check validity of command line parameters:
  if ( !isValid(static_cast<ETriangulationMethod>(triangulation_method))) {
    OPENMVG_LOG_ERROR << "Invalid triangulation method";
    return EXIT_FAILURE;
  }

  if ( !isValid(openMVG::cameras::EINTRINSIC(user_camera_model)) )  {
    OPENMVG_LOG_ERROR << "Invalid camera type";
    return EXIT_FAILURE;
  }

  const cameras::Intrinsic_Parameter_Type intrinsic_refinement_options =
      cameras::StringTo_Intrinsic_Parameter_Type(sIntrinsic_refinement_options);
  if (intrinsic_refinement_options == static_cast<cameras::Intrinsic_Parameter_Type>(0) )
  {
    OPENMVG_LOG_ERROR << "Invalid input for Bundle Adjustment Intrinsic parameter refinement option";
    return EXIT_FAILURE;
  }

  const sfm::Extrinsic_Parameter_Type extrinsic_refinement_options =
      sfm::StringTo_Extrinsic_Parameter_Type(sExtrinsic_refinement_options);
  if (extrinsic_refinement_options == static_cast<sfm::Extrinsic_Parameter_Type>(0) )
  {
    OPENMVG_LOG_ERROR << "Invalid input for the Bundle Adjustment Extrinsic parameter refinement option";
    return EXIT_FAILURE;
  }

  ESfMSceneInitializer scene_initializer_enum;
  if (!StringToEnum(sfm_initializer_method, scene_initializer_enum))
  {
    OPENMVG_LOG_ERROR << "Invalid input for the SfM initializer option";
    return EXIT_FAILURE;
  }

  ESfMEngine sfm_engine_type;
  if (!StringToEnum(engine_name, sfm_engine_type))
  {
    OPENMVG_LOG_ERROR << "Invalid input for the SfM Engine type";
    return EXIT_FAILURE;
  }

  if (rotation_averaging_method < ROTATION_AVERAGING_L1 ||
      rotation_averaging_method > ROTATION_AVERAGING_L2 )  {
    OPENMVG_LOG_ERROR << "Rotation averaging method is invalid";
    return EXIT_FAILURE;
  }

  // #ifndef USE_PATENTED_LIGT
  //   if (translation_averaging_method == TRANSLATION_LIGT) {
  //     OPENMVG_LOG_ERROR << "OpenMVG was not compiled with USE_PATENTED_LIGT cmake option";
  //     return EXIT_FAILURE;
  //   }
  // #endif
  //   if (translation_averaging_method < TRANSLATION_AVERAGING_L1 ||
  //       translation_averaging_method > TRANSLATION_LIGT )  {
  //     OPENMVG_LOG_ERROR << "Translation averaging method is invalid";
  //     return EXIT_FAILURE;
  //   }

  if (directory_output.empty())  {
    OPENMVG_LOG_ERROR << "It is an invalid output directory";
    return EXIT_FAILURE;
  }

  // SfM related

  // Load input SfM_Data scene
  SfM_Data sfm_data;
  const ESfM_Data sfm_data_loading_etypes =
      scene_initializer_enum == ESfMSceneInitializer::INITIALIZE_EXISTING_POSES ?
        ESfM_Data(VIEWS|INTRINSICS|EXTRINSICS) : ESfM_Data(VIEWS|INTRINSICS);
  if (!Load(sfm_data, filename_sfm_data, sfm_data_loading_etypes)) {
    OPENMVG_LOG_ERROR << "The input SfM_Data file \""<< filename_sfm_data << "\" cannot be read.";
    return EXIT_FAILURE;
  }

  if (!stlplus::folder_exists(directory_output))
  {
    if (!stlplus::folder_create(directory_output))
    {
      OPENMVG_LOG_ERROR << "Cannot create the output directory";
      return EXIT_FAILURE;
    }
  }

  //
  // Match and features
  //
  if (directory_match.empty() && !filename_match.empty() && stlplus::file_exists(filename_match))
  {
    directory_match = stlplus::folder_part(filename_match);
    filename_match = stlplus::filename_part(filename_match);
  }

  // Init the regions_type from the image describer file (used for image regions extraction)
  using namespace openMVG::features;
  const std::string sImage_describer = stlplus::create_filespec(directory_match, "image_describer", "json");
  std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);
  if (!regions_type)
  {
    OPENMVG_LOG_ERROR << "Invalid: " << sImage_describer << " regions type file.";
    return EXIT_FAILURE;
  }

  // Features reading
  std::shared_ptr<Features_Provider> feats_provider = std::make_shared<Features_Provider>();
  if (!feats_provider->load(sfm_data, directory_match, regions_type)) {
    OPENMVG_LOG_ERROR << "Cannot load view corresponding features in directory: " << directory_match << ".";
    return EXIT_FAILURE;
  }
  // Matches reading
  std::shared_ptr<Matches_Provider> matches_provider = std::make_shared<Matches_Provider>();
  if // Try to read the provided match filename or the default one (matches.f.txt/bin)
  (
  !(matches_provider->load(sfm_data, stlplus::create_filespec(directory_match, filename_match)) ||
      matches_provider->load(sfm_data, stlplus::create_filespec(directory_match, "matches.f.txt")) ||
      matches_provider->load(sfm_data, stlplus::create_filespec(directory_match, "matches.f.bin")) ||
      matches_provider->load(sfm_data, stlplus::create_filespec(directory_match, "matches.e.txt")) ||
      matches_provider->load(sfm_data, stlplus::create_filespec(directory_match, "matches.e.bin")))
      )
  {
    OPENMVG_LOG_ERROR << "Cannot load the match file.";
    return EXIT_FAILURE;
  }

  std::unique_ptr<SfMSceneInitializer> scene_initializer;
  switch(scene_initializer_enum)
  {
  case ESfMSceneInitializer::INITIALIZE_AUTO_PAIR:
    OPENMVG_LOG_ERROR << "Not yet implemented.";
    return EXIT_FAILURE;
    break;
  case ESfMSceneInitializer::INITIALIZE_MAX_PAIR:
    scene_initializer.reset(new SfMSceneInitializerMaxPair(sfm_data,
                                 feats_provider.get(),
                                 matches_provider.get()));
    break;
  case ESfMSceneInitializer::INITIALIZE_EXISTING_POSES:
    scene_initializer.reset(new SfMSceneInitializer(sfm_data,
                            feats_provider.get(),
                            matches_provider.get()));
    break;
  case ESfMSceneInitializer::INITIALIZE_STELLAR:
    scene_initializer.reset(new SfMSceneInitializerStellar(sfm_data,
                                 feats_provider.get(),
                                 matches_provider.get()));
    break;
  default:
    OPENMVG_LOG_ERROR << "Unknown SFM Scene initializer method";
    return EXIT_FAILURE;
  }
  if (!scene_initializer)
  {
    OPENMVG_LOG_ERROR << "Invalid scene initializer.";
    return EXIT_FAILURE;
  }


  std::unique_ptr<ReconstructionEngine> sfm_engine;
  switch (sfm_engine_type)
  {
  case ESfMEngine::INCREMENTAL:
  {
    SequentialSfMReconstructionEngine * engine =
        new SequentialSfMReconstructionEngine(
          sfm_data,
          directory_output,
          stlplus::create_filespec(directory_output, "Reconstruction_Report.html"));

    // Configuration:
    engine->SetFeaturesProvider(feats_provider.get());
    engine->SetMatchesProvider(matches_provider.get());

    // Configure reconstruction parameters
    engine->SetUnknownCameraType(EINTRINSIC(user_camera_model));
    engine->Set_Use_Motion_Prior(b_use_motion_priors);
    engine->SetTriangulationMethod(static_cast<ETriangulationMethod>(triangulation_method));
    engine->SetResectionMethod(static_cast<resection::SolverType>(resection_method));

    // Handle Initial pair parameter
    if (!initial_pair_string.first.empty() && !initial_pair_string.second.empty())
    {
      Pair initial_pair_index;
      if (!computeIndexFromImageNames(sfm_data, initial_pair_string, initial_pair_index))
      {
        OPENMVG_LOG_ERROR << "Could not find the initial pairs <" << initial_pair_string.first
                  <<  ", " << initial_pair_string.second << ">!";
        return EXIT_FAILURE;
      }
      engine->setInitialPair(initial_pair_index);
    }

    sfm_engine.reset(engine);
  }
    break;
  case ESfMEngine::INCREMENTALV2:
  {
    SequentialSfMReconstructionEngine2 * engine =
        new SequentialSfMReconstructionEngine2(
          scene_initializer.get(),
          sfm_data,
          directory_output,
          stlplus::create_filespec(directory_output, "Reconstruction_Report.html"));

    // Configuration:
    engine->SetFeaturesProvider(feats_provider.get());
    engine->SetMatchesProvider(matches_provider.get());

    // Configure reconstruction parameters
    engine->Set_Intrinsics_Refinement_Type(intrinsic_refinement_options);
    engine->SetUnknownCameraType(EINTRINSIC(user_camera_model));
    engine->Set_Use_Motion_Prior(b_use_motion_priors);
    engine->SetTriangulationMethod(static_cast<ETriangulationMethod>(triangulation_method));
    engine->SetResectionMethod(static_cast<resection::SolverType>(resection_method));

    sfm_engine.reset(engine);
  }
    break;
  case ESfMEngine::GLOBAL:
  {
    GlobalSfMReconstructionEngine_RelativeMotions * engine =
        new GlobalSfMReconstructionEngine_RelativeMotions(
          sfm_data,
          directory_output,
          stlplus::create_filespec(directory_output, "Reconstruction_Report.html"));

    // Configuration:
    engine->SetFeaturesProvider(feats_provider.get());
    engine->SetMatchesProvider(matches_provider.get());

    // Configure reconstruction parameters
    engine->Set_Intrinsics_Refinement_Type(intrinsic_refinement_options);
    engine->Set_Use_Motion_Prior(b_use_motion_priors);

    // Configure motion averaging method
    engine->SetRotationAveragingMethod(ERotationAveragingMethod(rotation_averaging_method));
    engine->SetTranslationAveragingMethod(ETranslationAveragingMethod(translation_averaging_method));

    sfm_engine.reset(engine);
  }
    break;
  default:
    break;
  }
  if (!sfm_engine)
  {
    OPENMVG_LOG_ERROR << "Cannot create the requested SfM Engine.";
    return EXIT_FAILURE;
  }

  sfm_engine->Set_Intrinsics_Refinement_Type(intrinsic_refinement_options);
  sfm_engine->Set_Extrinsics_Refinement_Type(extrinsic_refinement_options);

  //---------------------------------------
  // Sequential reconstruction process
  //---------------------------------------

  openMVG::system::Timer timer;

  if (sfm_engine->Process())
  {
    OPENMVG_LOG_INFO << " Total Sfm took (s): " << timer.elapsed();

    OPENMVG_LOG_INFO << "...Generating SfM_Report.html";
    Generate_SfM_Report(sfm_engine->Get_SfM_Data(),
              stlplus::create_filespec(directory_output, "SfMReconstruction_Report.html"));

    //-- Export to disk computed scene (data & viewable results)
    OPENMVG_LOG_INFO << "...Export SfM_Data to disk.";
    Save(sfm_engine->Get_SfM_Data(),
       stlplus::create_filespec(directory_output, "sfm_data", ".bin"),
       ESfM_Data(ALL));

    Save(sfm_engine->Get_SfM_Data(),
       stlplus::create_filespec(directory_output, "cloud_and_poses", ".ply"),
       ESfM_Data(ALL));

    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}

//conversion dep
  bool exportToOpenMVS(
    const SfM_Data & sfm_data,
    const std::string & sOutFile,
    const std::string & sOutDir,
    const int iNumThreads = 0
    )
  {
    // Create undistorted images directory structure
    if (!stlplus::is_folder(sOutDir))
    {
      stlplus::folder_create(sOutDir);
      if (!stlplus::is_folder(sOutDir))
      {
        OPENMVG_LOG_ERROR << "Cannot access to one of the desired output directory";
        return false;
      }
    }
    const std::string sOutSceneDir = stlplus::folder_part(sOutFile);
    const std::string sOutImagesDir = stlplus::folder_to_relative_path(sOutSceneDir, sOutDir);

    // Export data :
    _INTERFACE_NAMESPACE::Interface scene;
    size_t nPoses(0);
    const uint32_t nViews((uint32_t)sfm_data.GetViews().size());

    system::LoggerProgress my_progress_bar(nViews,"- PROCESS VIEWS -");

    // OpenMVG can have not contiguous index, use a map to create the required OpenMVS contiguous ID index
    std::map<openMVG::IndexT, uint32_t> map_intrinsic, map_view;

    // define a platform with all the intrinsic group
    for (const auto& intrinsic: sfm_data.GetIntrinsics())
    {
      if (isPinhole(intrinsic.second->getType()))
      {
        const Pinhole_Intrinsic * cam = dynamic_cast<const Pinhole_Intrinsic*>(intrinsic.second.get());
        if (map_intrinsic.count(intrinsic.first) == 0)
          map_intrinsic.insert(std::make_pair(intrinsic.first, scene.platforms.size()));
        _INTERFACE_NAMESPACE::Interface::Platform platform;
        // add the camera
        _INTERFACE_NAMESPACE::Interface::Platform::Camera camera;
        camera.width = cam->w();
        camera.height = cam->h();
        camera.K = cam->K();
        // sub-pose
        camera.R = Mat3::Identity();
        camera.C = Vec3::Zero();
        platform.cameras.push_back(camera);
        scene.platforms.push_back(platform);
      }
    }

    // define images & poses
    scene.images.reserve(nViews);

    const std::string mask_filename_global = stlplus::create_filespec(sfm_data.s_root_path, "mask", "png");

    for (const auto& view : sfm_data.GetViews())
    {
      ++my_progress_bar;
      const std::string srcImage = stlplus::create_filespec(sfm_data.s_root_path, view.second->s_Img_path);
      const std::string mask_filename_local = stlplus::create_filespec(sfm_data.s_root_path, stlplus::basename_part(srcImage) + "_mask", "png");
      const std::string maskName = stlplus::create_filespec(sOutImagesDir, stlplus::basename_part(srcImage) + ".mask.png");
      const std::string globalMaskName = stlplus::create_filespec(sOutImagesDir, "global_mask_" + std::to_string(view.second.get()->id_intrinsic), ".png");
      if (!stlplus::is_file(srcImage))
      {
        OPENMVG_LOG_INFO << "Cannot read the corresponding image: " << srcImage;
        return false;
      }

      if (sfm_data.IsPoseAndIntrinsicDefined(view.second.get())) 
      {
        map_view[view.first] = scene.images.size();

        _INTERFACE_NAMESPACE::Interface::Image image;
        image.name = stlplus::create_filespec(sOutImagesDir, view.second->s_Img_path);
        image.platformID = map_intrinsic.at(view.second->id_intrinsic);
        _INTERFACE_NAMESPACE::Interface::Platform& platform = scene.platforms[image.platformID];
        image.cameraID = 0;
        if (stlplus::file_exists(mask_filename_local))
          image.maskName = maskName;
        else if (stlplus::file_exists(mask_filename_global))
          image.maskName = globalMaskName;

        _INTERFACE_NAMESPACE::Interface::Platform::Pose pose;
        image.poseID = platform.poses.size();
        image.ID = map_view[view.first];
        const openMVG::geometry::Pose3 poseMVG(sfm_data.GetPoseOrDie(view.second.get()));
        pose.R = poseMVG.rotation();
        pose.C = poseMVG.center();
        platform.poses.push_back(pose);
        ++nPoses;

        scene.images.emplace_back(image);
      }
      else
      {
        OPENMVG_LOG_INFO << "Cannot read the corresponding pose or intrinsic of view " << view.first;
      }
    }

    // Export undistorted images
    system::LoggerProgress my_progress_bar_images(sfm_data.views.size(), "- UNDISTORT IMAGES " );
    std::atomic<bool> bOk(true); // Use a boolean to track the status of the loop process
  #ifdef OPENMVG_USE_OPENMP
    const unsigned int nb_max_thread = (iNumThreads > 0)? iNumThreads : omp_get_max_threads();

    #pragma omp parallel for schedule(dynamic) num_threads(nb_max_thread)
  #endif
    for (int i = 0; i < static_cast<int>(sfm_data.views.size()); ++i)
    {
      ++my_progress_bar_images;

      if (!bOk)
        continue;

      Views::const_iterator iterViews = sfm_data.views.begin();
      std::advance(iterViews, i);
      const View * view = iterViews->second.get();

      // Get image paths
      const std::string srcImage = stlplus::create_filespec(sfm_data.s_root_path, view->s_Img_path);
      const std::string imageName = stlplus::create_filespec(sOutDir, view->s_Img_path);
      const std::string mask_filename_local = stlplus::create_filespec(sfm_data.s_root_path, stlplus::basename_part(srcImage) + "_mask", "png");
      const std::string maskName = stlplus::create_filespec(sOutDir, stlplus::basename_part(srcImage) + ".mask.png");


      if (sfm_data.IsPoseAndIntrinsicDefined(view))
      {
        // export undistorted images
        const openMVG::cameras::IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view->id_intrinsic).get();
        if (cam->have_disto())
        {
          // undistort image and save it
          Image<openMVG::image::RGBColor> imageRGB, imageRGB_ud;
          Image<uint8_t> image_gray, image_gray_ud;
          try
          {
            if (ReadImage(srcImage.c_str(), &imageRGB))
            {
              UndistortImage(imageRGB, cam, imageRGB_ud, BLACK);
              bOk = bOk & WriteImage(imageName.c_str(), imageRGB_ud);
            }
            else // If RGBColor reading fails, try to read as gray image
            if (ReadImage(srcImage.c_str(), &image_gray))
            {
              UndistortImage(image_gray, cam, image_gray_ud, BLACK);
              const bool bRes = WriteImage(imageName.c_str(), image_gray_ud);
              bOk = bOk & bRes;
            }
            else
            {
              bOk = bOk & false;
            }

            Image<unsigned char> imageMask;
            // Try to read the local mask
            if (stlplus::file_exists(mask_filename_local))
            {
              if (!ReadImage(mask_filename_local.c_str(), &imageMask)||
                  !(imageMask.Width() == cam->w() && imageMask.Height() == cam->h()))
              {
                OPENMVG_LOG_ERROR
                  << "Invalid mask: " << mask_filename_local << ';';
                bOk = bOk & false;
                continue;
              }
              UndistortImage(imageMask, cam, image_gray_ud, BLACK);
              const bool bRes = WriteImage(maskName.c_str(), image_gray_ud);
              bOk = bOk & bRes;
            }
          }
          catch (const std::bad_alloc& e)
          {
            bOk = bOk & false;
          }
        }
        else
        {
          // just copy image
          stlplus::file_copy(srcImage, imageName);
          if (stlplus::file_exists(mask_filename_local))
          {
            stlplus::file_copy(mask_filename_local, maskName);
          }
        }
      }
      else
      {
        // just copy the image
        stlplus::file_copy(srcImage, imageName);
        if (stlplus::file_exists(mask_filename_local))
        {
            stlplus::file_copy(mask_filename_local, maskName);
        }
      }
    }
    if (stlplus::file_exists(mask_filename_global))
    {
      for (int i = 0; i < static_cast<int>(sfm_data.GetIntrinsics().size()); i++)
      {
        const openMVG::cameras::IntrinsicBase* cam = sfm_data.GetIntrinsics().at(i).get();
        const std::string maskName = stlplus::create_filespec(sOutDir, "global_mask_" + std::to_string(i), ".png");
        Image<uint8_t> imageMask;
        Image<uint8_t> image_gray, image_gray_ud;
        if (cam->have_disto())
        {
          // Try to read the global mask
          if (!ReadImage(mask_filename_global.c_str(), &imageMask)||
              !(imageMask.Width() == cam->w() && imageMask.Height() == cam->h()))
          {
            OPENMVG_LOG_ERROR
              << "Invalid global mask: " << mask_filename_global << ';';
            bOk = bOk & false;
          }
            UndistortImage(imageMask, cam, image_gray_ud, BLACK);
            const bool bRes = WriteImage(maskName.c_str(), image_gray_ud);
            bOk = bOk & bRes;
        }
        else
        {
          stlplus::file_copy(mask_filename_global, maskName);
        }
      }
    }

    if (!bOk)
    {
      OPENMVG_LOG_ERROR << "Catched a memory error in the image conversion."
      << " Please consider to use less threads ([-n|--numThreads])." << std::endl;
      return false;
    }

    // define structure
    scene.vertices.reserve(sfm_data.GetLandmarks().size());
    for (const auto& vertex: sfm_data.GetLandmarks())
    {
      const Landmark & landmark = vertex.second;
      _INTERFACE_NAMESPACE::Interface::Vertex vert;
      _INTERFACE_NAMESPACE::Interface::Vertex::ViewArr& views = vert.views;
      for (const auto& observation: landmark.obs)
      {
        const auto it(map_view.find(observation.first));
        if (it != map_view.end()) {
          _INTERFACE_NAMESPACE::Interface::Vertex::View view;
          view.imageID = it->second;
          view.confidence = 0;
          views.push_back(view);
        }
      }
      if (views.size() < 2)
        continue;
      std::sort(
        views.begin(), views.end(),
        [] (const _INTERFACE_NAMESPACE::Interface::Vertex::View& view0, const _INTERFACE_NAMESPACE::Interface::Vertex::View& view1)
        {
          return view0.imageID < view1.imageID;
        }
      );
      vert.X = landmark.X.cast<float>();
      scene.vertices.push_back(vert);
    }

    // write OpenMVS data
    if (!_INTERFACE_NAMESPACE::ARCHIVE::SerializeSave(scene, sOutFile))
      return false;

    OPENMVG_LOG_INFO
      << "Scene saved to OpenMVS interface format:\n"
      << " #platforms: " << scene.platforms.size();
      for (int i = 0; i < scene.platforms.size(); ++i)
      {
        OPENMVG_LOG_INFO << "  platform ( " << i << " ) #cameras: " << scene.platforms[i].cameras.size();
      }
    OPENMVG_LOG_INFO
      << "  " << scene.images.size() << " images (" << nPoses << " calibrated)\n"
      << "  " << scene.vertices.size() << " Landmarks";
    return true;
  }
//

int mvgToMvs(std::string sSfM_Data_Filename, std::string sOutFile, 
             std::string sOutDir)
{
  // CmdLine cmd;
  // std::string sSfM_Data_Filename;
  // std::string sOutFile = "scene.mvs";
  // std::string sOutDir = "undistorted_images";
  int iNumThreads = 0;

  // cmd.add( make_option('i', sSfM_Data_Filename, "sfmdata") );
  // cmd.add( make_option('o', sOutFile, "outfile") );
  // cmd.add( make_option('d', sOutDir, "outdir") );
  #ifdef OPENMVG_USE_OPENMP
    cmd.add( make_option('n', iNumThreads, "numThreads") );
  #endif

  // try {
  //   if (argc == 1) throw std::string("Invalid command line parameter.");
  //   cmd.process(argc, argv);
  // } catch (const std::string& s) {
  //   OPENMVG_LOG_INFO << "Usage: " << argv[0] << '\n'
  //     << "[-i|--sfmdata] filename, the SfM_Data file to convert\n"
  //     << "[-o|--outfile] OpenMVS scene file\n"
  //     << "[-d|--outdir] undistorted images path\n"
  // #ifdef OPENMVG_USE_OPENMP
  //       << "[-n|--numThreads] number of thread(s)\n"
  // #endif
  //     ;

  //   OPENMVG_LOG_ERROR << s;
  //   return EXIT_FAILURE;
  // }

  if (stlplus::extension_part(sOutFile) != "mvs") {
    OPENMVG_LOG_ERROR
      << "Invalid output file extension: " << sOutFile << "."
      << "You must use a filename with a .mvs extension.";
      return EXIT_FAILURE;
  }

  // Read the input SfM scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(ALL))) {
    OPENMVG_LOG_ERROR << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read.";
    return EXIT_FAILURE;
  }

  // Export OpenMVS data structure
  if (!exportToOpenMVS(sfm_data, sOutFile, sOutDir, iNumThreads))
  {
    OPENMVG_LOG_ERROR << "The output openMVS scene file cannot be written";
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}