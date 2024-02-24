#include <string>
int initImageList(std::string sImageDir, std::string sOutputDir, std::string sfileDatabase);
//1. Compute features                openMVG_main_ComputeFeatures  
int getFeatures(std::string sSfM_Data_Filename, std::string sOutDir);
//2. Compute pairs                   openMVG_main_PairGenerator
int imagePairer(std::string sSfMDataFilename, std::string sOutputPairsFilename,
                std::string sPairMode = "EXHAUSTIVE", int iContiguousCount = -1);
//3. Compute matches                 openMVG_main_ComputeMatches
//BREAKS UNDER OPENMP...
int imageMatcher(std::string  sSfM_Data_Filename, std::string sOutputMatchesFilename, std::string  sPredefinedPairList,
                 float fDistRatio = 0.8f, std::string sNearestMatchingMethod = "AUTO", 
                 bool bForce = false, unsigned int ui_max_cache_size = 0);
//4. Filter matches                  openMVG_main_GeometricFilter
int imageFilterGeo(std::string sSfM_Data_Filename, std::string sPutativeMatchesFilename, 
                    std::string sFilteredMatchesFilename);
//5. Incremental reconstruction      openMVG_main_SfM
int imageSfm(std::string filename_sfm_data, std::string directory_match, std::string directory_output,
             std::string engine_name = "INCREMENTAL");
//11. Export to openMVS              openMVG_main_openMVG2openMVS
int mvgToMvs(std::string sSfM_Data_Filename, std::string sOutFile = "scene.mvs", 
             std::string sOutDir = "undistorted_images");