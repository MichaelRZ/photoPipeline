#include "mvs.h"
// #undef _USE_OPENMP

//densify pointcloud
#include "MVS/Common.h"
#include "MVS/Scene.h"
#include <boost/program_options.hpp>
#include <vector>

using namespace MVS;
#define APPNAME _T("Densify")
//densify pointcloud
  namespace {

  namespace OPT {
  String strInputFileName;
  String strOutputFileName;
  String strViewNeighborsFileName;
  String strOutputViewNeighborsFileName;
  String strMeshFileName;
  String strExportROIFileName;
  String strImportROIFileName;
  String strDenseConfigFileName;
  String strExportDepthMapsName;
  float fMaxSubsceneArea;
  float fSampleMesh;
  float fBorderROI;
  bool bCrop2ROI;
  int nEstimateROI;
  int nFusionMode;
  int thFilterPointCloud;
  int nExportNumViews;
  int nArchiveType;
  int nProcessPriority;
  unsigned nMaxThreads;
  String strConfigFileName;
  boost::program_options::variables_map vm;
  } // namespace OPT

  // initialize and parse the command line parameters
  bool Initialize(size_t argc, LPCTSTR* argv)
  {
    // initialize log and console
    OPEN_LOG();
    OPEN_LOGCONSOLE();

    // group of options allowed only on command line
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,h", "produce this help message")
      ("working-folder,w", boost::program_options::value<std::string>(&WORKING_FOLDER), "working directory (default current directory)")
      ("config-file,c", boost::program_options::value<std::string>(&OPT::strConfigFileName)->default_value(APPNAME _T(".cfg")), "file name containing program options")
      ("archive-type", boost::program_options::value(&OPT::nArchiveType)->default_value(ARCHIVE_DEFAULT), "project archive type: 0-text, 1-binary, 2-compressed binary")
      ("process-priority", boost::program_options::value(&OPT::nProcessPriority)->default_value(-1), "process priority (below normal by default)")
      ("max-threads", boost::program_options::value(&OPT::nMaxThreads)->default_value(0), "maximum number of threads (0 for using all available cores)")
      #if TD_VERBOSE != TD_VERBOSE_OFF
      ("verbosity,v", boost::program_options::value(&g_nVerbosityLevel)->default_value(
        #if TD_VERBOSE == TD_VERBOSE_DEBUG
        3
        #else
        2
        #endif
        ), "verbosity level")
      #endif
      #ifdef _USE_CUDA
      ("cuda-device", boost::program_options::value(&CUDA::desiredDeviceID)->default_value(-1), "CUDA device number to be used for depth-map estimation (-2 - CPU processing, -1 - best GPU, >=0 - device index)")
      #endif
      ;

    // group of options allowed both on command line and in config file
    #ifdef _USE_CUDA
    const unsigned nNumViewsDefault(8);
    const unsigned numIters(4);
    #else
    const unsigned nNumViewsDefault(5);
    const unsigned numIters(3);
    #endif
    unsigned nResolutionLevel;
    unsigned nMaxResolution;
    unsigned nMinResolution;
    unsigned nNumViews;
    unsigned nMinViewsFuse;
    unsigned nSubResolutionLevels;
    unsigned nEstimationIters;
    unsigned nEstimationGeometricIters;
    unsigned nEstimateColors;
    unsigned nEstimateNormals;
    unsigned nOptimize;
    int nIgnoreMaskLabel;
    bool bRemoveDmaps;
    boost::program_options::options_description config("Densify options");
    config.add_options()
      ("input-file,i", boost::program_options::value<std::string>(&OPT::strInputFileName), "input filename containing camera poses and image list")
      ("output-file,o", boost::program_options::value<std::string>(&OPT::strOutputFileName), "output filename for storing the dense point-cloud (optional)")
      ("view-neighbors-file", boost::program_options::value<std::string>(&OPT::strViewNeighborsFileName), "input filename containing the list of views and their neighbors (optional)")
      ("output-view-neighbors-file", boost::program_options::value<std::string>(&OPT::strOutputViewNeighborsFileName), "output filename containing the generated list of views and their neighbors")
      ("resolution-level", boost::program_options::value(&nResolutionLevel)->default_value(1), "how many times to scale down the images before point cloud computation")
      ("max-resolution", boost::program_options::value(&nMaxResolution)->default_value(2560), "do not scale images higher than this resolution")
      ("min-resolution", boost::program_options::value(&nMinResolution)->default_value(640), "do not scale images lower than this resolution")
      ("sub-resolution-levels", boost::program_options::value(&nSubResolutionLevels)->default_value(2), "number of patch-match sub-resolution iterations (0 - disabled)")
      ("number-views", boost::program_options::value(&nNumViews)->default_value(nNumViewsDefault), "number of views used for depth-map estimation (0 - all neighbor views available)")
      ("number-views-fuse", boost::program_options::value(&nMinViewsFuse)->default_value(3), "minimum number of images that agrees with an estimate during fusion in order to consider it inlier (<2 - only merge depth-maps)")
      ("ignore-mask-label", boost::program_options::value(&nIgnoreMaskLabel)->default_value(-1), "integer value for the label to ignore in the segmentation mask (<0 - disabled)")
      ("iters", boost::program_options::value(&nEstimationIters)->default_value(numIters), "number of patch-match iterations")
      ("geometric-iters", boost::program_options::value(&nEstimationGeometricIters)->default_value(2), "number of geometric consistent patch-match iterations (0 - disabled)")
      ("estimate-colors", boost::program_options::value(&nEstimateColors)->default_value(2), "estimate the colors for the dense point-cloud (0 - disabled, 1 - final, 2 - estimate)")
      ("estimate-normals", boost::program_options::value(&nEstimateNormals)->default_value(2), "estimate the normals for the dense point-cloud (0 - disabled, 1 - final, 2 - estimate)")
      ("sub-scene-area", boost::program_options::value(&OPT::fMaxSubsceneArea)->default_value(0.f), "split the scene in sub-scenes such that each sub-scene surface does not exceed the given maximum sampling area (0 - disabled)")
      ("sample-mesh", boost::program_options::value(&OPT::fSampleMesh)->default_value(0.f), "uniformly samples points on a mesh (0 - disabled, <0 - number of points, >0 - sample density per square unit)")
      ("fusion-mode", boost::program_options::value(&OPT::nFusionMode)->default_value(0), "depth-maps fusion mode (-2 - fuse disparity-maps, -1 - export disparity-maps only, 0 - depth-maps & fusion, 1 - export depth-maps only)")
      ("postprocess-dmaps", boost::program_options::value(&nOptimize)->default_value(7), "flags used to filter the depth-maps after estimation (0 - disabled, 1 - remove-speckles, 2 - fill-gaps, 4 - adjust-filter)")
      ("filter-point-cloud", boost::program_options::value(&OPT::thFilterPointCloud)->default_value(0), "filter dense point-cloud based on visibility (0 - disabled)")
      ("export-number-views", boost::program_options::value(&OPT::nExportNumViews)->default_value(0), "export points with >= number of views (0 - disabled, <0 - save MVS project too)")
          ("roi-border", boost::program_options::value(&OPT::fBorderROI)->default_value(0), "add a border to the region-of-interest when cropping the scene (0 - disabled, >0 - percentage, <0 - absolute)")
          ("estimate-roi", boost::program_options::value(&OPT::nEstimateROI)->default_value(2), "estimate and set region-of-interest (0 - disabled, 1 - enabled, 2 - adaptive)")
          ("crop-to-roi", boost::program_options::value(&OPT::bCrop2ROI)->default_value(true), "crop scene using the region-of-interest")
          ("remove-dmaps", boost::program_options::value(&bRemoveDmaps)->default_value(false), "remove depth-maps after fusion")
          ;

    // hidden options, allowed both on command line and
    // in config file, but will not be shown to the user
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("mesh-file", boost::program_options::value<std::string>(&OPT::strMeshFileName), "mesh file name used for image pair overlap estimation")
      ("export-roi-file", boost::program_options::value<std::string>(&OPT::strExportROIFileName), "ROI file name to be exported form the scene")
      ("import-roi-file", boost::program_options::value<std::string>(&OPT::strImportROIFileName), "ROI file name to be imported into the scene")
      ("dense-config-file", boost::program_options::value<std::string>(&OPT::strDenseConfigFileName), "optional configuration file for the densifier (overwritten by the command line options)")
      ("export-depth-maps-name", boost::program_options::value<std::string>(&OPT::strExportDepthMapsName), "render given mesh and save the depth-map for every image to this file name base (empty - disabled)")
      ;

    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(config).add(hidden);

    boost::program_options::options_description config_file_options;
    config_file_options.add(config).add(hidden);
    std::cout << "before input" << std::endl;
    boost::program_options::positional_options_description p;
    p.add("input-file", -1);

    try {
      // parse command line options
      boost::program_options::store(boost::program_options::command_line_parser((int)argc, argv).options(cmdline_options).positional(p).run(), OPT::vm);
      boost::program_options::notify(OPT::vm);
      INIT_WORKING_FOLDER;
      // parse configuration file
      std::ifstream ifs(MAKE_PATH_SAFE(OPT::strConfigFileName));
      if (ifs) {
        boost::program_options::store(parse_config_file(ifs, config_file_options), OPT::vm);
        boost::program_options::notify(OPT::vm);
      }
    }
    catch (const std::exception& e) {
      LOG(e.what());
      return false;
    }

    std::cout << "log file" << std::endl;

    // initialize the log file
    OPEN_LOGFILE(MAKE_PATH(APPNAME _T("-")+Util::getUniqueName(0)+_T(".log")).c_str());

    // print application details: version and command line
    Util::LogBuild();
    LOG(_T("Command line: ") APPNAME _T("%s"), Util::CommandLineToString(argc, argv).c_str());

    std::cout << "vlidte input" << std::endl;

    // validate input
    Util::ensureValidPath(OPT::strInputFileName);
    if (OPT::vm.count("help") || OPT::strInputFileName.empty()) {
      boost::program_options::options_description visible("Available options");
      visible.add(generic).add(config);
      GET_LOG() << visible;
    }
    if (OPT::strInputFileName.empty())
      return false;

    // initialize optional options
    Util::ensureValidPath(OPT::strOutputFileName);
    Util::ensureValidPath(OPT::strViewNeighborsFileName);
    Util::ensureValidPath(OPT::strOutputViewNeighborsFileName);
    Util::ensureValidPath(OPT::strMeshFileName);
    Util::ensureValidPath(OPT::strExportROIFileName);
    Util::ensureValidPath(OPT::strImportROIFileName);
    if (OPT::strOutputFileName.empty())
      OPT::strOutputFileName = Util::getFileFullName(OPT::strInputFileName) + _T("_dense.mvs");

    // init dense options
    if (!OPT::strDenseConfigFileName.empty())
      OPT::strDenseConfigFileName = MAKE_PATH_SAFE(OPT::strDenseConfigFileName);
    OPTDENSE::init();
    const bool bValidConfig(OPTDENSE::oConfig.Load(OPT::strDenseConfigFileName));
    OPTDENSE::update();
    OPTDENSE::nResolutionLevel = nResolutionLevel;
    OPTDENSE::nMaxResolution = nMaxResolution;
    OPTDENSE::nMinResolution = nMinResolution;
    OPTDENSE::nSubResolutionLevels = nSubResolutionLevels;
    OPTDENSE::nNumViews = nNumViews;
    OPTDENSE::nMinViewsFuse = nMinViewsFuse;
    OPTDENSE::nEstimationIters = nEstimationIters;
    OPTDENSE::nEstimationGeometricIters = nEstimationGeometricIters;
    OPTDENSE::nEstimateColors = nEstimateColors;
    OPTDENSE::nEstimateNormals = nEstimateNormals;
    OPTDENSE::nOptimize = nOptimize;
    OPTDENSE::nIgnoreMaskLabel = nIgnoreMaskLabel;
    OPTDENSE::bRemoveDmaps = bRemoveDmaps;
    if (!bValidConfig && !OPT::strDenseConfigFileName.empty())
      OPTDENSE::oConfig.Save(OPT::strDenseConfigFileName);

    std::cout << "before openmp/brkpd" << std::endl;

    // initialize global options
    Process::setCurrentProcessPriority((Process::Priority)OPT::nProcessPriority);
    #ifdef _USE_OPENMP
    if (OPT::nMaxThreads != 0)
      omp_set_num_threads(OPT::nMaxThreads);
    #endif

    #ifdef _USE_BREAKPAD
    // start memory dumper
    MiniDumper::Create(APPNAME, WORKING_FOLDER);
    #endif

    Util::Init();
    std::cout << "end init" << std::endl;

    return true;
  }

  // finalize application instance
  void Finalize()
  {
    #if TD_VERBOSE != TD_VERBOSE_OFF
    // print memory statistics
    Util::LogMemoryInfo();
    #endif

    CLOSE_LOGFILE();
    CLOSE_LOGCONSOLE();
    CLOSE_LOG();
  }

  } // unnamed namespace

  int densePtCld(std::vector<std::string> arguments)
  {
    #ifdef _DEBUGINFO
      std::cout<< "debug info" << std::endl;
    #endif
    #ifdef _USE_CUDA
      std::cout << "use cuda" << std::endl;
    #endif
    #ifdef _USE_BREAKPAD
      std::cout << "use breakpad" << std::endl;
    #endif
    #ifdef _USE_OPENMP
      std::cout << "use openmp" << std::endl;
    #endif

    // std::vector<std::string> arguments = {"--dir", "/some_path"};
    std::vector<const char*> argv;
    for (const auto& arg : arguments)
        argv.push_back((const char*)arg.data());
    argv.push_back(nullptr);

    #ifdef _DEBUGINFO
    // set _crtBreakAlloc index to stop in <dbgheap.c> at allocation
    _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);// | _CRTDBG_CHECK_ALWAYS_DF);
    #endif

    // if (!Initialize(argc, argv))
    //   return EXIT_FAILURE;
    if (!Initialize(argv.size() - 1, argv.data()))
      return EXIT_FAILURE;
    std::cout << "threads" << std::endl;
    std::cout << OPT::nMaxThreads << std::endl;

    MVS::Scene scene(OPT::nMaxThreads);
    if (OPT::fSampleMesh != 0) {
      // sample input mesh and export the obtained point-cloud
      if (!scene.Load(MAKE_PATH_SAFE(OPT::strInputFileName), true) || scene.mesh.IsEmpty())
        return EXIT_FAILURE;
      TD_TIMER_START();
      PointCloud pointcloud;
      if (OPT::fSampleMesh > 0)
        scene.mesh.SamplePoints(OPT::fSampleMesh, 0, pointcloud);
      else
        scene.mesh.SamplePoints(ROUND2INT<unsigned>(-OPT::fSampleMesh), pointcloud);
      VERBOSE("Sample mesh completed: %u points (%s)", pointcloud.GetSize(), TD_TIMER_GET_FMT().c_str());
      pointcloud.Save(MAKE_PATH_SAFE(Util::getFileFullName(OPT::strOutputFileName))+_T(".ply"));
      Finalize();
      return EXIT_SUCCESS;
    }
    std::cout << "past first conditional" << std::endl;

    // load and estimate a dense point-cloud
    if (!scene.Load(MAKE_PATH_SAFE(OPT::strInputFileName)))
      return EXIT_FAILURE;
    std::cout << "past load" << std::endl;
    if (!OPT::strImportROIFileName.empty()) {
      std::ifstream fs(MAKE_PATH_SAFE(OPT::strImportROIFileName));
      if (!fs)
        return EXIT_FAILURE;
      fs >> scene.obb;
      scene.Save(MAKE_PATH_SAFE(Util::getFileFullName(OPT::strOutputFileName))+_T(".mvs"), (ARCHIVE_TYPE)OPT::nArchiveType);
      Finalize();
      return EXIT_SUCCESS;
    }
    std::cout << "past if" << std::endl;
    if (!scene.IsBounded())
      scene.EstimateROI(OPT::nEstimateROI, 1.1f);
    std::cout << "if 2" << std::endl;
    if (!OPT::strExportROIFileName.empty() && scene.IsBounded()) {
      std::ofstream fs(MAKE_PATH_SAFE(OPT::strExportROIFileName));
      if (!fs)
        return EXIT_FAILURE;
      fs << scene.obb;
      Finalize();
      return EXIT_SUCCESS;
    }
    std::cout << "if 3" << std::endl;
    if (!OPT::strMeshFileName.empty())
      scene.mesh.Load(MAKE_PATH_SAFE(OPT::strMeshFileName));
    if (!OPT::strViewNeighborsFileName.empty())
      scene.LoadViewNeighbors(MAKE_PATH_SAFE(OPT::strViewNeighborsFileName));
    if (!OPT::strOutputViewNeighborsFileName.empty()) {
      if (!scene.ImagesHaveNeighbors()) {
        VERBOSE("error: neighbor views not computed yet");
        return EXIT_FAILURE;
      }
      scene.SaveViewNeighbors(MAKE_PATH_SAFE(OPT::strOutputViewNeighborsFileName));
      return EXIT_SUCCESS;
    }
    if (!OPT::strExportDepthMapsName.empty() && !scene.mesh.IsEmpty()) {
      // project mesh onto each image and save the resulted depth-maps
      TD_TIMER_START();
      if (!scene.ExportMeshToDepthMaps(MAKE_PATH_SAFE(OPT::strExportDepthMapsName)))
        return EXIT_FAILURE;
      VERBOSE("Mesh projection completed: %u depth-maps (%s)", scene.images.size(), TD_TIMER_GET_FMT().c_str());
      Finalize();
      return EXIT_SUCCESS;
    }
    if (OPT::fMaxSubsceneArea > 0) {
      // split the scene in sub-scenes by maximum sampling area
      Scene::ImagesChunkArr chunks;
      scene.Split(chunks, OPT::fMaxSubsceneArea);
      scene.ExportChunks(chunks, GET_PATH_FULL(OPT::strOutputFileName), (ARCHIVE_TYPE)OPT::nArchiveType);
      Finalize();
      return EXIT_SUCCESS;
    }
    if (OPT::thFilterPointCloud < 0) {
      // filter point-cloud based on camera-point visibility intersections
      scene.PointCloudFilter(OPT::thFilterPointCloud);
      const String baseFileName(MAKE_PATH_SAFE(Util::getFileFullName(OPT::strOutputFileName))+_T("_filtered"));
      scene.Save(baseFileName+_T(".mvs"), (ARCHIVE_TYPE)OPT::nArchiveType);
      scene.pointcloud.Save(baseFileName+_T(".ply"));
      Finalize();
      return EXIT_SUCCESS;
    }
    if (OPT::nExportNumViews && scene.pointcloud.IsValid()) {
      // export point-cloud containing only points with N+ views
      const String baseFileName(MAKE_PATH_SAFE(Util::getFileFullName(OPT::strOutputFileName))+
        String::FormatString(_T("_%dviews"), ABS(OPT::nExportNumViews)));
      if (OPT::nExportNumViews > 0) {
        // export point-cloud containing only points with N+ views
        scene.pointcloud.SaveNViews(baseFileName+_T(".ply"), (IIndex)OPT::nExportNumViews);
      } else {
        // save scene and export point-cloud containing only points with N+ views
        scene.pointcloud.RemoveMinViews((IIndex)-OPT::nExportNumViews);
        scene.Save(baseFileName+_T(".mvs"), (ARCHIVE_TYPE)OPT::nArchiveType);
        scene.pointcloud.Save(baseFileName+_T(".ply"));
      }
      Finalize();
      return EXIT_SUCCESS;
    }
    if ((ARCHIVE_TYPE)OPT::nArchiveType != ARCHIVE_MVS) {
      #if TD_VERBOSE != TD_VERBOSE_OFF
      if (VERBOSITY_LEVEL > 1 && !scene.pointcloud.IsEmpty())
        scene.pointcloud.PrintStatistics(scene.images.data(), &scene.obb);
      #endif
      TD_TIMER_START();
      if (!scene.DenseReconstruction(OPT::nFusionMode, OPT::bCrop2ROI, OPT::fBorderROI)) {
        if (ABS(OPT::nFusionMode) != 1)
          return EXIT_FAILURE;
        VERBOSE("Depth-maps estimated (%s)", TD_TIMER_GET_FMT().c_str());
        Finalize();
        return EXIT_SUCCESS;
      }
      VERBOSE("Densifying point-cloud completed: %u points (%s)", scene.pointcloud.GetSize(), TD_TIMER_GET_FMT().c_str());
    }
    std::cout << "before saving" << std::endl;

    // save the final point-cloud
    const String baseFileName(MAKE_PATH_SAFE(Util::getFileFullName(OPT::strOutputFileName)));
    scene.Save(baseFileName+_T(".mvs"), (ARCHIVE_TYPE)OPT::nArchiveType);
    scene.pointcloud.Save(baseFileName+_T(".ply"));
    #if TD_VERBOSE != TD_VERBOSE_OFF
    if (VERBOSITY_LEVEL > 2)
      scene.ExportCamerasMLP(baseFileName+_T(".mlp"), baseFileName+_T(".ply"));
    #endif

    Finalize();
    return EXIT_SUCCESS;
  }
  /*----------------------------------------------------------------*/

//reconstruct mesh part

/*
 * ReconstructMesh.cpp
 *
 * Copyright (c) 2014-2015 SEACAVE
 *
 * Author(s):
 *
 *      cDc <cdc.seacave@gmail.com>
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Additional Terms:
 *
 *      You are required to preserve legal notices and author attributions in
 *      that material or in the Appropriate Legal Notices displayed by works
 *      containing it.
 */



// D E F I N E S ///////////////////////////////////////////////////
#undef APPNAME
#define APPNAME _T("ReconstructMesh")

// uncomment to enable multi-threading based on OpenMP
#ifdef _USE_OPENMP
#define RECMESH_USE_OPENMP
#endif


// S T R U C T S ///////////////////////////////////////////////////

namespace {

namespace OPT_mesh {
String strInputFileName;
String strOutputFileName;
String strMeshFileName;
bool bMeshExport;
float fDistInsert;
bool bUseOnlyROI;
bool bUseConstantWeight;
bool bUseFreeSpaceSupport;
float fThicknessFactor;
float fQualityFactor;
float fDecimateMesh;
unsigned nTargetFaceNum;
float fRemoveSpurious;
bool bRemoveSpikes;
unsigned nCloseHoles;
unsigned nSmoothMesh;
float fEdgeLength;
bool bCrop2ROI;
float fBorderROI;
float fSplitMaxArea;
unsigned nArchiveType;
int nProcessPriority;
unsigned nMaxThreads;
String strImagePointsFileName;
String strExportType;
String strConfigFileName;
boost::program_options::variables_map vm;
} // namespace OPT

// initialize and parse the command line parameters
bool Initialize_mesh(size_t argc, LPCTSTR* argv)
{
	// initialize log and console
	OPEN_LOG();
	OPEN_LOGCONSOLE();

	// group of options allowed only on command line
	boost::program_options::options_description generic("Generic options");
	generic.add_options()
		("help,h", "produce this help message")
		("working-folder,w", boost::program_options::value<std::string>(&WORKING_FOLDER), "working directory (default current directory)")
		("config-file,c", boost::program_options::value<std::string>(&OPT_mesh::strConfigFileName)->default_value(APPNAME _T(".cfg")), "file name containing program options")
		("export-type", boost::program_options::value<std::string>(&OPT_mesh::strExportType)->default_value(_T("ply")), "file type used to export the 3D scene (ply or obj)")
		("archive-type", boost::program_options::value(&OPT_mesh::nArchiveType)->default_value(ARCHIVE_DEFAULT), "project archive type: 0-text, 1-binary, 2-compressed binary")
		("process-priority", boost::program_options::value(&OPT_mesh::nProcessPriority)->default_value(-1), "process priority (below normal by default)")
		("max-threads", boost::program_options::value(&OPT_mesh::nMaxThreads)->default_value(0), "maximum number of threads (0 for using all available cores)")
		#if TD_VERBOSE != TD_VERBOSE_OFF
		("verbosity,v", boost::program_options::value(&g_nVerbosityLevel)->default_value(
			#if TD_VERBOSE == TD_VERBOSE_DEBUG
			3
			#else
			2
			#endif
			), "verbosity level")
		#endif
		#ifdef _USE_CUDA
		("cuda-device", boost::program_options::value(&CUDA::desiredDeviceID)->default_value(-1), "CUDA device number to be used to reconstruct the mesh (-2 - CPU processing, -1 - best GPU, >=0 - device index)")
		#endif
		;

	// group of options allowed both on command line and in config file
	boost::program_options::options_description config_main("Reconstruct options");
	config_main.add_options()
		("input-file,i", boost::program_options::value<std::string>(&OPT_mesh::strInputFileName), "input filename containing camera poses and image list")
		("output-file,o", boost::program_options::value<std::string>(&OPT_mesh::strOutputFileName), "output filename for storing the mesh")
		("min-point-distance,d", boost::program_options::value(&OPT_mesh::fDistInsert)->default_value(2.5f), "minimum distance in pixels between the projection of two 3D points to consider them different while triangulating (0 - disabled)")
		("integrate-only-roi", boost::program_options::value(&OPT_mesh::bUseOnlyROI)->default_value(false), "use only the points inside the ROI")
		("constant-weight", boost::program_options::value(&OPT_mesh::bUseConstantWeight)->default_value(true), "considers all view weights 1 instead of the available weight")
		("free-space-support,f", boost::program_options::value(&OPT_mesh::bUseFreeSpaceSupport)->default_value(false), "exploits the free-space support in order to reconstruct weakly-represented surfaces")
		("thickness-factor", boost::program_options::value(&OPT_mesh::fThicknessFactor)->default_value(1.f), "multiplier adjusting the minimum thickness considered during visibility weighting")
		("quality-factor", boost::program_options::value(&OPT_mesh::fQualityFactor)->default_value(1.f), "multiplier adjusting the quality weight considered during graph-cut")
		;
	boost::program_options::options_description config_clean("Clean options");
	config_clean.add_options()
		("decimate", boost::program_options::value(&OPT_mesh::fDecimateMesh)->default_value(1.f), "decimation factor in range (0..1] to be applied to the reconstructed surface (1 - disabled)")
		("target-face-num", boost::program_options::value(&OPT_mesh::nTargetFaceNum)->default_value(0), "target number of faces to be applied to the reconstructed surface. (0 - disabled)")
		("remove-spurious", boost::program_options::value(&OPT_mesh::fRemoveSpurious)->default_value(20.f), "spurious factor for removing faces with too long edges or isolated components (0 - disabled)")
		("remove-spikes", boost::program_options::value(&OPT_mesh::bRemoveSpikes)->default_value(true), "flag controlling the removal of spike faces")
		("close-holes", boost::program_options::value(&OPT_mesh::nCloseHoles)->default_value(30), "try to close small holes in the reconstructed surface (0 - disabled)")
		("smooth", boost::program_options::value(&OPT_mesh::nSmoothMesh)->default_value(2), "number of iterations to smooth the reconstructed surface (0 - disabled)")
		("edge-length", boost::program_options::value(&OPT_mesh::fEdgeLength)->default_value(0.f), "remesh such that the average edge length is this size (0 - disabled)")
		("roi-border", boost::program_options::value(&OPT_mesh::fBorderROI)->default_value(0), "add a border to the region-of-interest when cropping the scene (0 - disabled, >0 - percentage, <0 - absolute)")
		("crop-to-roi", boost::program_options::value(&OPT_mesh::bCrop2ROI)->default_value(true), "crop scene using the region-of-interest")
		;

	// hidden options, allowed both on command line and
	// in config file, but will not be shown to the user
	boost::program_options::options_description hidden("Hidden options");
	hidden.add_options()
		("mesh-file", boost::program_options::value<std::string>(&OPT_mesh::strMeshFileName), "mesh file name to clean (skips the reconstruction step)")
		("mesh-export", boost::program_options::value(&OPT_mesh::bMeshExport)->default_value(false), "just export the mesh contained in loaded project")
		("split-max-area", boost::program_options::value(&OPT_mesh::fSplitMaxArea)->default_value(0.f), "maximum surface area that a sub-mesh can contain (0 - disabled)")
		("image-points-file", boost::program_options::value<std::string>(&OPT_mesh::strImagePointsFileName), "input filename containing the list of points from an image to project on the mesh (optional)")
		;

	boost::program_options::options_description cmdline_options;
	cmdline_options.add(generic).add(config_main).add(config_clean).add(hidden);

	boost::program_options::options_description config_file_options;
	config_file_options.add(config_main).add(config_clean).add(hidden);

	boost::program_options::positional_options_description p;
	p.add("input-file", -1);

	try {
		// parse command line options
		boost::program_options::store(boost::program_options::command_line_parser((int)argc, argv).options(cmdline_options).positional(p).run(), OPT_mesh::vm);
		boost::program_options::notify(OPT_mesh::vm);
		INIT_WORKING_FOLDER;
		// parse configuration file
		std::ifstream ifs(MAKE_PATH_SAFE(OPT_mesh::strConfigFileName));
		if (ifs) {
			boost::program_options::store(parse_config_file(ifs, config_file_options), OPT_mesh::vm);
			boost::program_options::notify(OPT_mesh::vm);
		}
	}
	catch (const std::exception& e) {
		LOG(e.what());
		return false;
	}

	// initialize the log file
	OPEN_LOGFILE(MAKE_PATH(APPNAME _T("-")+Util::getUniqueName(0)+_T(".log")).c_str());

	// print application details: version and command line
	Util::LogBuild();
	LOG(_T("Command line: ") APPNAME _T("%s"), Util::CommandLineToString(argc, argv).c_str());

	// validate input
	Util::ensureValidPath(OPT_mesh::strInputFileName);
	Util::ensureUnifySlash(OPT_mesh::strInputFileName);
	if (OPT_mesh::vm.count("help") || OPT_mesh::strInputFileName.IsEmpty()) {
		boost::program_options::options_description visible("Available options");
		visible.add(generic).add(config_main).add(config_clean);
		GET_LOG() << visible;
	}
	if (OPT_mesh::strInputFileName.IsEmpty())
		return false;
	OPT_mesh::strExportType = OPT_mesh::strExportType.ToLower() == _T("obj") ? _T(".obj") : _T(".ply");

	// initialize optional options
	Util::ensureValidPath(OPT_mesh::strOutputFileName);
	Util::ensureUnifySlash(OPT_mesh::strOutputFileName);
	Util::ensureValidPath(OPT_mesh::strImagePointsFileName);
	if (OPT_mesh::strOutputFileName.IsEmpty())
		OPT_mesh::strOutputFileName = Util::getFileFullName(OPT_mesh::strInputFileName) + _T("_mesh.mvs");

	// initialize global options
	Process::setCurrentProcessPriority((Process::Priority)OPT_mesh::nProcessPriority);
	#ifdef _USE_OPENMP
	if (OPT_mesh::nMaxThreads != 0)
		omp_set_num_threads(OPT_mesh::nMaxThreads);
	#endif

	#ifdef _USE_BREAKPAD
	// start memory dumper
	MiniDumper::Create(APPNAME, WORKING_FOLDER);
	#endif

	Util::Init();
	return true;
}

// finalize application instance
void Finalize_mesh()
{
	#if TD_VERBOSE != TD_VERBOSE_OFF
	// print memory statistics
	Util::LogMemoryInfo();
	#endif

	CLOSE_LOGFILE();
	CLOSE_LOGCONSOLE();
	CLOSE_LOG();
}

} // unnamed namespace


// export 3D coordinates corresponding to 2D coordinates provided by inputFileName:
// parse image point list; first line is the name of the image to project,
// each consequent line store the xy coordinates to project:
// <image-name> <number-of-points>
// <x-coord1> <y-coord1>
// <x-coord2> <y-coord2>
// ...
// 
// for example:
// N01.JPG 3
// 3090 2680
// 3600 2100
// 3640 2190
bool Export3DProjections(Scene& scene, const String& inputFileName) {
	SML smlPointList(_T("ImagePoints"));
	smlPointList.Load(inputFileName);
	ASSERT(smlPointList.GetArrChildren().size() <= 1);
	IDX idx(0);

	// read image name
	size_t argc;
	CAutoPtrArr<LPSTR> argv;
	while (true) {
		argv = Util::CommandLineToArgvA(smlPointList.GetValue(idx).val, argc);
		if (argc > 0 && argv[0][0] != _T('#'))
			break;
		if (++idx == smlPointList.size())
			return false;
	}
	if (argc < 2)
		return false;
	String imgName(argv[0]);
	IIndex imgID(NO_ID);
	for (const Image& imageData : scene.images) {
		if (!imageData.IsValid())
			continue;
		if (imageData.name.substr(imageData.name.size() - imgName.size()) == imgName) {
			imgID = imageData.ID;
			break;
		}
	}
	if (imgID == NO_ID) {
		VERBOSE("Unable to find image named: %s", imgName.c_str());
		return false;
	}

	// read image points
	std::vector<Point2f> imagePoints;
	while (++idx != smlPointList.size()) {
		// parse image element
		const String& line(smlPointList.GetValue(idx).val);
		argv = Util::CommandLineToArgvA(line, argc);
		if (argc > 0 && argv[0][0] == _T('#'))
			continue;
		if (argc < 2) {
			VERBOSE("Invalid image coordinates: %s", line.c_str());
			continue;
		}
		const Point2f pt(
			String::FromString<float>(argv[0], -1),
			String::FromString<float>(argv[1], -1));
		if (pt.x > 0 && pt.y > 0)
			imagePoints.emplace_back(pt);
	}
	if (imagePoints.empty()) {
		VERBOSE("Unable to read image points from: %s", imgName.c_str());
		return false;
	}

	// prepare output file
	String outFileName(Util::insertBeforeFileExt(inputFileName, "_3D"));
	File oStream(outFileName, File::WRITE, File::CREATE | File::TRUNCATE);
	if (!oStream.isOpen()) {
		VERBOSE("Unable to open output file: %s", outFileName.c_str());
		return false;
	}

	// print image name
	oStream.print("%s %u\n", imgName.c_str(), imagePoints.size());

	// init mesh octree
	const Mesh::Octree octree(scene.mesh.vertices, [](Mesh::Octree::IDX_TYPE size, Mesh::Octree::Type /*radius*/) {
		return size > 256;
	});
	scene.mesh.ListIncidenteFaces();

	// save 3D coord in the output file
	const Image& imgToExport = scene.images[imgID];
	for (const Point2f& pt : imagePoints) {
		// define ray from camera center to each x,y image coord
		const Ray3 ray(imgToExport.camera.C, normalized(imgToExport.camera.RayPoint<REAL>(pt)));
		// find ray intersection with the mesh
		const IntersectRayMesh intRay(octree, ray, scene.mesh);
		if (intRay.pick.IsValid()) {
			const Point3d ptHit(ray.GetPoint(intRay.pick.dist));
			oStream.print("%.7f %.7f %.7f\n", ptHit.x, ptHit.y, ptHit.z);
		} else 
			oStream.print("NA\n");
	}
	return true;
}

int reconstructMesh(std::vector<std::string> arguments)
{
  std::vector<const char*> argv;
  for (const auto& arg : arguments)
      argv.push_back((const char*)arg.data());
  argv.push_back(nullptr);

	#ifdef _DEBUGINFO
	// set _crtBreakAlloc index to stop in <dbgheap.c> at allocation
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);// | _CRTDBG_CHECK_ALWAYS_DF);
	#endif
  std::cout << "mesh init" << std::endl;
	if (!Initialize_mesh(argv.size() - 1, argv.data()))
		return EXIT_FAILURE;
  std::cout << "scene" << std::endl;
	Scene scene(OPT_mesh::nMaxThreads);
  std::cout << "scene load" << std::endl;
	// load project
	if (!scene.Load(MAKE_PATH_SAFE(OPT_mesh::strInputFileName), OPT_mesh::fSplitMaxArea > 0 || OPT_mesh::fDecimateMesh < 1 || OPT_mesh::nTargetFaceNum > 0))
		return EXIT_FAILURE;
	const String baseFileName(MAKE_PATH_SAFE(Util::getFileFullName(OPT_mesh::strOutputFileName)));
	if (OPT_mesh::fSplitMaxArea > 0) {
		// split mesh using max-area constraint
		Mesh::FacesChunkArr chunks;
		if (scene.mesh.Split(chunks, OPT_mesh::fSplitMaxArea))
			scene.mesh.Save(chunks, baseFileName);
		Finalize_mesh();
		return EXIT_SUCCESS;
	}
  std::cout << "past max area" << std::endl;
	if (!OPT_mesh::strImagePointsFileName.empty() && !scene.mesh.IsEmpty()) {
		Export3DProjections(scene, MAKE_PATH_SAFE(OPT_mesh::strImagePointsFileName));
		return EXIT_SUCCESS;
	}
  std::cout << "past export 3d proj" << std::endl;
	if (OPT_mesh::bMeshExport) {
    std::cout << "if case" << std::endl;
		// check there is a mesh to export
		if (scene.mesh.IsEmpty())
			return EXIT_FAILURE;
		// save mesh
		const String fileName(MAKE_PATH_SAFE(OPT_mesh::strOutputFileName));
		scene.mesh.Save(fileName);
		#if TD_VERBOSE != TD_VERBOSE_OFF
		if (VERBOSITY_LEVEL > 2)
			scene.ExportCamerasMLP(baseFileName+_T(".mlp"), fileName);
		#endif
	} else {
    std::cout << "else case" << std::endl;

		const OBB3f initialOBB(scene.obb);
		if (OPT_mesh::fBorderROI > 0)
			scene.obb.EnlargePercent(OPT_mesh::fBorderROI);
		else if (OPT_mesh::fBorderROI < 0)
			scene.obb.Enlarge(-OPT_mesh::fBorderROI);
		if (OPT_mesh::strMeshFileName.IsEmpty() && scene.mesh.IsEmpty()) {
			// reset image resolution to the original size and
			// make sure the image neighbors are initialized before deleting the point-cloud
			#ifdef RECMESH_USE_OPENMP
			bool bAbort(false);
			#pragma omp parallel for
			for (int_t idx=0; idx<(int_t)scene.images.GetSize(); ++idx) {
				#pragma omp flush (bAbort)
				if (bAbort)
					continue;
				const uint32_t idxImage((uint32_t)idx);
			#else
			FOREACH(idxImage, scene.images) {
			#endif
				Image& imageData = scene.images[idxImage];
				if (!imageData.IsValid())
					continue;
				// reset image resolution
				if (!imageData.ReloadImage(0, false)) {
					#ifdef RECMESH_USE_OPENMP
					bAbort = true;
					#pragma omp flush (bAbort)
					continue;
					#else
					return EXIT_FAILURE;
					#endif
				}
				imageData.UpdateCamera(scene.platforms);
				// select neighbor views
				if (imageData.neighbors.IsEmpty()) {
					IndexArr points;
					scene.SelectNeighborViews(idxImage, points);
				}
			}
			#ifdef RECMESH_USE_OPENMP
			if (bAbort)
				return EXIT_FAILURE;
			#endif
			// reconstruct a coarse mesh from the given point-cloud
      std::cout << "here" << std::endl;
			TD_TIMER_START();
			if (OPT_mesh::bUseConstantWeight)
				scene.pointcloud.pointWeights.Release();
      std::cout << "crash??" << std::endl;
      //seems to crash around here......
			if (!scene.ReconstructMesh(OPT_mesh::fDistInsert, OPT_mesh::bUseFreeSpaceSupport, OPT_mesh::bUseOnlyROI, 4, OPT_mesh::fThicknessFactor, OPT_mesh::fQualityFactor))
				return EXIT_FAILURE;
      std::cout << "no crash" << std::endl;
			VERBOSE("Mesh reconstruction completed: %u vertices, %u faces (%s)", scene.mesh.vertices.GetSize(), scene.mesh.faces.GetSize(), TD_TIMER_GET_FMT().c_str());
			#if TD_VERBOSE != TD_VERBOSE_OFF
			if (VERBOSITY_LEVEL > 2) {
				// dump raw mesh
				scene.mesh.Save(baseFileName+_T("_raw")+OPT_mesh::strExportType);
			}
			#endif
		} else if (!OPT_mesh::strMeshFileName.IsEmpty()) {
			// load existing mesh to clean
			scene.mesh.Load(MAKE_PATH_SAFE(OPT_mesh::strMeshFileName));
		}
    std::cout << "clean mesh" << std::endl;
		// clean the mesh
		if (OPT_mesh::bCrop2ROI && scene.IsBounded()) {
			TD_TIMER_START();
			const size_t numVertices = scene.mesh.vertices.size();
			const size_t numFaces = scene.mesh.faces.size();
			scene.mesh.RemoveFacesOutside(scene.obb);
			VERBOSE("Mesh trimmed to ROI: %u vertices and %u faces removed (%s)",
				numVertices-scene.mesh.vertices.size(), numFaces-scene.mesh.faces.size(), TD_TIMER_GET_FMT().c_str());
		}
		const float fDecimate(OPT_mesh::nTargetFaceNum ? static_cast<float>(OPT_mesh::nTargetFaceNum) / scene.mesh.faces.size() : OPT_mesh::fDecimateMesh);
		scene.mesh.Clean(fDecimate, OPT_mesh::fRemoveSpurious, OPT_mesh::bRemoveSpikes, OPT_mesh::nCloseHoles, OPT_mesh::nSmoothMesh, OPT_mesh::fEdgeLength, false);
		scene.mesh.Clean(1.f, 0.f, OPT_mesh::bRemoveSpikes, OPT_mesh::nCloseHoles, 0u, 0.f, false); // extra cleaning trying to close more holes
		scene.mesh.Clean(1.f, 0.f, false, 0u, 0u, 0.f, true); // extra cleaning to remove non-manifold problems created by closing holes
		scene.obb = initialOBB;

    std::cout << "save mesh" << std::endl;
		// save the final mesh
		scene.Save(baseFileName+_T(".mvs"), (ARCHIVE_TYPE)OPT_mesh::nArchiveType);
		scene.mesh.Save(baseFileName+OPT_mesh::strExportType);
		#if TD_VERBOSE != TD_VERBOSE_OFF
		if (VERBOSITY_LEVEL > 2)
			scene.ExportCamerasMLP(baseFileName+_T(".mlp"), baseFileName+OPT_mesh::strExportType);
		#endif
	}
  std::cout << "save 1" << std::endl;
	if (!OPT_mesh::strImagePointsFileName.empty()) {
		Export3DProjections(scene, MAKE_PATH_SAFE(OPT_mesh::strImagePointsFileName));
		return EXIT_SUCCESS;
	}
  std::cout << "before finalize" << std::endl;
	Finalize_mesh();
	return EXIT_SUCCESS;
}
/*----------------------------------------------------------------*/
