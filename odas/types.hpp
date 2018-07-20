#ifndef __ODAS_TYPES
#define __ODAS_TYPES

#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/octree/octree.h>

#define PointType pcl::PointXYZ
#define PointOCT pcl::PointXYZL
#define PCDT pcl::PointCloud<pcl::PointXYZ>
#define PCDLT pcl::PointCloud<pcl::PointXYZL>
#define KdtreeT pcl::KdTreeFLANN<pcl::PointXYZ>
#define OctreeT pcl::octree::OctreePointCloudPointVector<pcl::PointXYZ>

#define LIDAR_ANGLE (10.0f / 180.0f * M_PI)

#define _USE_MATH_DEFINES
//#define TEST_VERBOSE
#define TEST_BENCHMARK
//#define TEST_SAVE_OCT

#endif
