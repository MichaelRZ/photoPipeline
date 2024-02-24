cmake_minimum_required(VERSION 3.10)

set (CMAKE_CXX_STANDARD 14)

set(CMAKE_TOOLCHAIN_FILE ../vcpkg/scripts/buildsystems/vcpkg.cmake)
# set(CMAKE_CXX_FLAGS_RELEASE "-O1")

project(imitation)
# list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_LIST_DIR}/../vcpkg/installed/x64-osx/")
# list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/../vcpkg/installed/x64-osx/share/")

# MACRO(register_definitions DEF)
#   add_definitions(${DEF})
#   string(REPLACE "-D" "" STRIPPED_DEF ${DEF})
#   list(APPEND OpenMVG_DEFINITIONS ${STRIPPED_DEF})
# ENDMACRO()
# register_definitions(-DOPENMVG_STD_UNORDERED_MAP)

# find_package(OpenMVG CONFIG REQUIRED PATHS /Users/michaelr/Desktop/vcpkg/installed/x64-osx/ NO_DEFAULT_PATH)

find_package(OpenMVG CONFIG REQUIRED)
add_library(mvg mvg.cpp)
target_include_directories(mvg PRIVATE ${OPENMVG_INCLUDE_DIRS})
target_include_directories(mvg PRIVATE ../nonFree/sift/vl)
target_compile_definitions(mvg PRIVATE -DOPENMVG_STD_UNORDERED_MAP)
target_compile_options(mvg PRIVATE -fvisibility=hidden)
#https://stackoverflow.com/questions/20643368/cmake-different-include-directories-for-different-targets

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1 -mavx2")
# note: 16 additional targets are not displayed.
  # target_link_libraries(main PRIVATE OpenMVG::openMVG_sfm OpenMVG::openMVG_exif OpenMVG::openMVG_fast OpenMVG::openMVG_kvld)

# include_directories(${VLFEAT_INCLUDE_DIRS})

# include_directories(${PROJECT_SOURCE_DIR})
# SET(GCC_COVERAGE_COMPILE_FLAGS "-g")
# add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})
# find_package(cereal CONFIG REQUIRED)
add_executable(main main.cpp)
# target_compile_options(main PRIVATE -fvisibility=default )

# add_executable(main main.cpp)
add_subdirectory(nonFree)
# find_package(vlfeat CONFIG REQUIRED)

#possible openmp
# if(NOT TARGET OpenMP::OpenMP_CXX)
#     find_package(Threads REQUIRED)
#     add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
#     set_property(TARGET OpenMP::OpenMP_CXX
#                  PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
#     # Only works if the same flag is passed to the linker; use CMake 3.9+ otherwise (Intel, AppleClang)
#     set_property(TARGET OpenMP::OpenMP_CXX
#                  PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)

# endif()

#actual openmp
# find_package(OpenMP REQUIRED)
# add_compile_definitions(OPENMVG_USE_OPENMP)
# target_link_libraries(main PUBLIC OpenMP::OpenMP_CXX)

#openmvs
# this is heuristically generated, and may not be correct
# list(APPEND CMAKE_PREFIX_PATH "Users/michaelr/Desktop/imitation/vcg/")
# find_package(zstd CONFIG REQUIRED)
# target_link_libraries(main PRIVATE $<IF:$<TARGET_EXISTS:zstd::libzstd_shared>,zstd::libzstd_shared,zstd::libzstd_static>)
# https://stackoverflow.com/questions/70710800/how-to-connect-mysql-to-c-on-m1-mac/70907595#70907595
#export LIBRARY_PATH=$LIBRARY_PATH:$(brew --prefix zstd)/lib/
set(CMAKE_CXX_VISIBILITY_PRESET default)
set(CMAKE_C_VISIBILITY_PRESET default)
# https://www.alecjacobson.com/weblog/2578.html
# find_package(COLMAP REQUIRED)
# target_link_libraries(main colmap::colmap)
# find_package(OpenMVS CONFIG REQUIRED PATHS /Users/michaelr/Desktop/vcpkg/installed/x64-osx)
# find_package(Eigen3 REQUIRED PATHS /Users/michaelr/Desktop/eigen-3.2.10)
find_package(OpenMVS CONFIG REQUIRED)
add_library(mvs mvs.cpp)
target_compile_options(mvs PRIVATE -fvisibility=default)

target_include_directories(mvs PRIVATE ${OpenMVS_INCLUDE_DIRS})

# include_directories(../colmap/src/)
include(CMakePrintHelpers)
cmake_print_variables(OpenMVS_DEFINITIONS)
cmake_print_variables(OpenMVS_INCLUDE_DIRS)
target_compile_definitions(mvs PRIVATE ${OpenMVS_DEFINITIONS})
# add_definitions(${OpenMVS_DEFINITIONS})
# target_link_libraries(main PRIVATE ${OpenMVS_LIBRAIRES})
# target_link_libraries(main PRIVATE 
#     # OpenMVS::IO 
#     OpenMVS::MVS 
#     # OpenMVS::Math 
#     # OpenMVS::Common
#     # colmap::colmap
#     # colmap::colmap_controllers
#     # colmap::colmap_estimators
#     # colmap::colmap_exe
#     # colmap::colmap_feature
#     # colmap::colmap_geometry
#     # colmap::colmap_image
#     # colmap::colmap_math
#     # colmap::colmap_mvs
#     # colmap::colmap_optim
#     # colmap::colmap_retrieval
#     # colmap::colmap_scene
#     # colmap::colmap_sensor
#     # colmap::colmap_sfm
#     # colmap::colmap_util
#     # # Third-party.
#     # colmap::colmap_lsd
#     # colmap::colmap_poisson_recon
#     # colmap::colmap_vlfeat
# )
target_link_libraries(mvg
  PRIVATE
    OpenMVG::openMVG_camera
    OpenMVG::openMVG_exif
    OpenMVG::openMVG_geodesy
    OpenMVG::openMVG_image
    OpenMVG::openMVG_numeric
    OpenMVG::openMVG_sfm
    OpenMVG::openMVG_system
    OpenMVG::openMVG_features
    OpenMVG::openMVG_matching
    OpenMVG::openMVG_matching_image_collection
    OpenMVG::openMVG_geometry

    # OpenMVG::openMVG_camera
    # OpenMVG::openMVG_exif
    # OpenMVG::openMVG_features
    # OpenMVG::openMVG_geodesy
    # OpenMVG::openMVG_geometry
    # OpenMVG::openMVG_graph
    # OpenMVG::openMVG_image
    # OpenMVG::openMVG_linearProgramming
    # OpenMVG::openMVG_matching
    # OpenMVG::openMVG_matching_image_collection
    # OpenMVG::openMVG_multiview
    # OpenMVG::openMVG_numeric
    # OpenMVG::openMVG_robust_estimation
    # OpenMVG::openMVG_sfm
    # OpenMVG::openMVG_system
    # OpenMVS::MVS
    vlsift
)
target_link_libraries(mvs PRIVATE OpenMVS::MVS)
target_link_libraries(main PRIVATE mvs mvg)
# target_link_libraries(main PRIVATE mvs)
