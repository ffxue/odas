#ifndef __ODAS_SYMPCD
#define __ODAS_SYMPCD


#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/common/common.h>
#include <pcl/point_types.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/octree/octree.h>
#include <pcl/filters/voxel_grid.h>
//#include <pcl/keypoints/sift_keypoint.h>
//#include <pcl/features/normal_3d.h>

#include "types.hpp"

class SymPCD
{
private:
    std::string filename;
    PCDT::Ptr   cloud;
    PCDT::Ptr   oct_voxel_cloud;
    PCDLT::Ptr  oct_search;
    KdtreeT     kdtree, octkdtree;
    // boundary and distances
    PointType   pmin, pmax, pcenter;
    float       dist_diag, dist_diag_xy, p_dist_avg, dist_apr, dist_half_edge;
    float       sq_dist_apr, sq_dist_diag, epsilon;
    // z-slices
    std::vector<PCDT::Ptr> z_slices;
    std::vector<KdtreeT> z_kdtrees;
    float       z_slice_dense, z_slice_sparse, z_roof_est;
    int         base_z_dense;
    // memory for kdtree search
    std::vector<int> idx_search;
    std::vector<float> squared_dist;
    // sift
    PCDT::Ptr   sift_cloud;
    KdtreeT     sift_kdtree;

    float M_PI_HALF;
    float M_PI_ONE_AND_HALF;

    void prepare_boundary()
    {
        // boundary
        pcl::getMinMax3D (*cloud, pmin, pmax);
        pcenter.x = (pmin.x + pmax.x) / 2;
        pcenter.y = (pmin.y + pmax.y) / 2;
        pcenter.z = (pmin.z + pmax.z) / 2;
        
        // set up distances
        float width = pmax.x - pmin.x;
        float depth = pmax.y - pmin.y;
        float height = pmax.z - pmin.z;
        sq_dist_diag = width * width + depth * depth + height * height;
        dist_diag = std::sqrt(sq_dist_diag);
        dist_diag_xy = std::sqrt(width * width + depth * depth);
        //dist_half_edge = std::min(std::min(width, depth), height) / 2;
        dist_half_edge = std::min(width, depth) / 2;
        // effective_area = top face (1) + projected side faces (4)
        float effective_area = width * depth + 2 * (width + depth) * height 
                                     * std::tan(LIDAR_ANGLE);
        p_dist_avg = std::sqrt(effective_area / cloud->size());
        dist_apr = dist_diag * epsilon;
        sq_dist_apr = dist_apr * dist_apr;
        
#ifdef TEST_VERBOSE
        std::cerr << ">> density = " << cloud->size() / effective_area << ", svg_dist = " << p_dist_avg << ", dist_apr = " << dist_apr << std::endl;
#endif
    }
    
    void global_shift_xy()
    {
        if (cloud->size() == 0)
            return;
        
        for (auto it = cloud->begin(); it != cloud->end(); ++it)
        {
            it->x -= pcenter.x;
            it->y -= pcenter.y;
        }
    }

    void prepare_octree(int depth)
    {
        // filtering for search pcd by octree
        OctreeT octree(dist_half_edge / std::pow(2, depth-1));
        //OctreeT octree(dist_diag_xy / std::pow(2, depth));
        octree.setInputCloud(cloud);
        octree.defineBoundingBox();
        octree.addPointsFromInputCloud();
        oct_voxel_cloud->clear();
        oct_voxel_cloud->reserve(octree.getLeafCount());
        
        Eigen::Vector3f min_pt, max_pt;
        for (auto it = octree.leaf_begin(); it != octree.leaf_end(); it++)
        {
            // get voxel center
            octree.getVoxelBounds(it, min_pt, max_pt);
            PointType p;
            p.x = (min_pt[0]+max_pt[0])/2;
            p.y = (min_pt[1]+max_pt[1])/2;
            p.z = (min_pt[2]+max_pt[2])/2;
            // snap to real cloud
            if (kdtree.nearestKSearch(p, 1, idx_search, squared_dist) > 0)
            {
                p.x = (*cloud)[idx_search[0]].x;
                p.y = (*cloud)[idx_search[0]].y;
                p.z = (*cloud)[idx_search[0]].z;
            }
            // form a XYZ cloud
            oct_voxel_cloud->push_back(p);
            // form a XYZL cloud, L = weight (number of representing points)
            PointOCT p1;
            p1.x = p.x;
            p1.y = p.y;
            p1.z = p.z;
            p1.label = it.getLeafContainer().getSize();
            oct_search->push_back(p1);
        }
        octree.deleteTree();
        octkdtree.setInputCloud(oct_voxel_cloud);
    }

    void prepare_z_slices()
    {
        // set up shared variables
        z_slice_dense = p_dist_avg;
        z_slice_sparse = p_dist_avg / std::tan(LIDAR_ANGLE);
        
        int height_bins = (int)std::floor((pmax.z - pmin.z) / p_dist_avg)+1;
        int* height_sum = (int*)std::malloc(height_bins * sizeof(int));
        std::memset(height_sum, 0, height_bins);
        for (auto it = cloud->begin(); it != cloud->end(); ++it)
        {
            int target_bin = (int)std::floor((it->z - pmin.z) / p_dist_avg);
            height_sum[target_bin]++;
        }
        int num_roof_p_est = (int) ((pmax.x - pmin.x) * (pmax.y - pmin.y) 
                                / (p_dist_avg * p_dist_avg));
        int accumulated_p = 0;
        for (int i = height_bins - 1; i >= 0; i--)
        {
            z_roof_est = pmin.z + p_dist_avg * i;
            accumulated_p += height_sum[i];
            if (accumulated_p > num_roof_p_est)
                break;
        }
        // free mem
        free(height_sum);
        
        // do z-slicing
        base_z_dense = static_cast<int>(floor((z_roof_est - pmin.z) / z_slice_sparse)) 
                        + 1;
        int num_bins =  base_z_dense
                        + static_cast<int>(floor((pmax.z - z_roof_est) / z_slice_dense))
                        + 1;
        for (int i = 0; i < num_bins; i++)
        {
            PCDT::Ptr pc(new PCDT);
            z_slices.push_back(pc);
        }
        // setup slices
        for (auto it = cloud->begin(); it != cloud->end(); ++it)
        {
            PointType p(*it);
            z_slices[zidx(it->z)]->push_back(p);
        }
        // extract kd-trees
        for (int i = 0; i < num_bins; i++)
        {
            if (z_slices[i]->size() > 0)
            {
                //pcl::PointCloud<PointType>::Ptr pc1(new pcl::PointCloud<PointType>);
                //pcl::VoxelGrid<PointType> sor;
                //sor.setInputCloud (slices[i]);
                //sor.setLeafSize (p_dist_avg, p_dist_avg, p_dist_avg);
                //sor.filter(*pc1);

                KdtreeT tree;
                tree.setInputCloud(z_slices[i]);
                z_kdtrees.push_back(tree);
            }
            else
            {
                z_kdtrees.push_back(NULL);
            }
        }
    }
    
    void prepare_sift()
    {
        /*
        // Parameters for sift computation
        const float min_scale = std::max(0.01f, p_dist_avg);
        const int n_octaves = 3;
        const int n_scales_per_octave = 4;
        const float min_contrast = 0.001f;

        // Estimate the normals of the cloud
        pcl::NormalEstimation<PointType, pcl::PointNormal> ne;
        pcl::PointCloud<pcl::PointNormal>::Ptr cloud_normals (new pcl::PointCloud<pcl::PointNormal>);
        pcl::search::KdTree<PointType>::Ptr tree_n(new pcl::search::KdTree<PointType>());

        ne.setInputCloud(cloud);
        ne.setSearchMethod(tree_n);
        ne.setRadiusSearch(min_scale * 5);
        ne.compute(*cloud_normals);

        // Copy the xyz info from cloud and add it to cloud_normals as the xyz field in PointNormals estimation is zero
        for(size_t i = 0; i<cloud_normals->points.size(); ++i)
        {
            cloud_normals->points[i].x = cloud->points[i].x;
            cloud_normals->points[i].y = cloud->points[i].y;
            cloud_normals->points[i].z = cloud->points[i].z;
        }

        // Estimate the sift interest points using normals values from xyz as the Intensity variants
        pcl::SIFTKeypoint<pcl::PointNormal, pcl::PointWithScale> sift;
        pcl::PointCloud<pcl::PointWithScale> result;
        pcl::search::KdTree<pcl::PointNormal>::Ptr tree(new pcl::search::KdTree<pcl::PointNormal> ());
        sift.setSearchMethod(tree);
        sift.setScales(min_scale, n_octaves, n_scales_per_octave);
        sift.setMinimumContrast(min_contrast);
        sift.setInputCloud(cloud_normals);
        sift.compute(result);
        pcl::copyPointCloud(result, *sift_cloud);
        sift_kdtree.setInputCloud(sift_cloud);
        */
    }


public:
    
    SymPCD()
        : cloud(new PCDT), oct_voxel_cloud(new PCDT), oct_search(new PCDLT)
            , sift_cloud(new PCDT), idx_search(1), squared_dist(1)
    {
        M_PI_HALF = M_PI/2.0f;
        M_PI_ONE_AND_HALF = M_PI + M_PI/2.0f;
    }
    
    void load(std::string fname, int oct_depth = 4, float epsilon_diag = 0.005)
    {
        // load file
        filename = fname;
        std::string fileext = fname.substr(fname.size() - 3);
        if (fileext == "pcd")
        {
            pcl::PCDReader reader;
            reader.read (filename.c_str(), *cloud);
        }
        else if (fileext == "ply")
        {
            pcl::PLYReader reader;
            reader.read (filename.c_str(), *cloud);
        }
        else
        {
            std::cerr << "[Error] Unknown extension \"." << fileext 
                    << "\". (.pcd or .ply expected)" << std::endl;
            exit(0);
        }


        epsilon = epsilon_diag;
        prepare_boundary();
        global_shift_xy();
        
        kdtree.setInputCloud (cloud);
        prepare_octree(oct_depth);
        prepare_z_slices();
#ifdef TEST_SAVE_OCT
        std::stringstream ss;
        ss << "octree" << filename << "_" << oct_depth << ".pcd";
        pcl::io::savePCDFileASCII (ss.str(), *oct_voxel_cloud);
#endif
    }

    void get_reflection(PointType *pout, const float x, const float y, const float rho, const float theta)
    {
        if (theta == M_PI_HALF)
        {
            pout->x = x;
            pout->y = rho + rho - y;
        }
        else if (theta == M_PI_ONE_AND_HALF)
        {
            pout->x = x;
            pout->y = - rho - rho - y;
        }
        else
        {
            float dx = 2 * (rho - std::cos(theta) * x - std::sin(theta) * y ) * std::cos(theta);
            pout->x = x + dx;
            pout->y = y + dx * std::tan(theta);
        }
        
    }

    double fast_approx_fitness(const float rho, const float theta)
    {
        float squared_error = 0.0f;
        int sum_apr = 0;
        bool apr;
        PointType p;
        auto octpc_end = oct_search->end();
        for (auto it = oct_search->begin(); it != octpc_end; it++)
        {
            // error metrics
            get_reflection(&p, it->x, it->y, rho, theta);
            p.z = it->z;
            squared_error += get_approx_mixed_penalty(p, &apr) * it->label;  // times weight
            if (apr)
                sum_apr += it->label;
            //squared_error += get_approx_nn_sq_dist(p);
        }
        return (double)std::sqrt(squared_error / cloud->size())     // rmse
                + dist_diag * (cloud->size() - sum_apr) / cloud->size();
        //return 1.0 - (double)found_apr*1.0 / cloud->size();
    }
    
    int zidx(float z)
    {
        // using dense slices for roofs
        if (z >= z_roof_est)
            return static_cast<int>(floor((z - z_roof_est) / z_slice_dense)) + base_z_dense;
        return static_cast<int>(floor((z - pmin.z) / z_slice_sparse));
    }

    float get_approx_mixed_penalty(const PointType p, bool* apr_ptr)
    {
        *apr_ptr = false;
        if (z_kdtrees[zidx(p.z)].nearestKSearch(p, 1, idx_search, squared_dist) > 0)
        {
            if (squared_dist[0] < sq_dist_apr)
                *apr_ptr = true;
            return squared_dist[0];
        }
        return sq_dist_diag;
    }

    float get_approx_nn_sq_dist(const PointType p)
    {
        //int zidx = static_cast<int>(floor((p.z - pmin.z) / z_slice_sparse));
        if (octkdtree.nearestKSearch(p, 1, idx_search, squared_dist) > 0)
        {
            return squared_dist[0];
        }
        return sq_dist_diag;
    }

    int get_approx_apr(const PointType p)
    {
        int zidx = static_cast<int>(floor((p.z - pmin.z) / z_slice_sparse));
        if (z_kdtrees[zidx].nearestKSearch(p, 1, idx_search, squared_dist) > 0)
            if (std::sqrt(squared_dist[0]) < dist_apr)
                return 1;
        return 0;
    }

    double validate(const float ro, const float theta, double * apr = nullptr)
    {
        int   all_cnt = 0;
        float all_sum = 0.0f;
        float all_squ = 0.0f;

        PointType p;
        auto end_of_cloud = cloud->end();
        for (auto it = cloud->begin(); it != end_of_cloud; ++it)
        {
            get_reflection(&p, it->x, it->y, ro, theta);
            p.z = it->z;
            if (kdtree.nearestKSearch(p, 1, idx_search, squared_dist) > 0)
            {
                all_sum += std::sqrt(squared_dist[0]);
                all_squ += squared_dist[0];
                if (squared_dist[0] < sq_dist_apr)
                    all_cnt++;
            }
        }
#ifdef TEST_VERBOSE
        std::cerr << "EVALUATION: apr=" << all_cnt*1.0/cloud->size() << ", rmse=" << std::sqrt(all_squ*1.0/cloud->size()) << std::endl;
#endif
        if (apr != nullptr)
            *apr = all_cnt*1.0/cloud->size();
        return (double)std::sqrt(all_squ / cloud->size());
    }
    
    double segment(const float ro, const float theta)
    {
        PCDT::Ptr   pc_left(new PCDT), pc_right(new PCDT), pc_rem(new PCDT);

        PointType p;
        auto end_of_cloud = cloud->end();
        for (auto it = cloud->begin(); it != end_of_cloud; ++it)
        {
            get_reflection(&p, it->x, it->y, ro, theta);
            p.z = it->z;
            if (kdtree.nearestKSearch(p, 1, idx_search, squared_dist) > 0)
            {
                PointType p0(it->x, it->y, it->z);
                if (squared_dist[0] < sq_dist_apr)
                {
                    // corresponded
                    if (p.x >= it->x)
                    {
                        // left cloud
                        pc_left->push_back(p0);
                    }
                    else
                    {
                        // right cloud
                        pc_right->push_back(p0);
                    }
                }
                else
                {
                    // add to rem
                    pc_rem->push_back(p0);
                }
            }
        }
        
        std::stringstream ss;
        if (pc_left->size() > 0)
        {
            ss << "segments/" << filename << "_left.pcd";
            std::cout << ss.str() << std::endl;
            pcl::io::savePCDFileBinary (ss.str(), *pc_left);
        }
        if (pc_right->size() > 0)
        {
            ss.str("");
            ss << "segments/" << filename << "_right.pcd";
            std::cout << ss.str() << std::endl;
            pcl::io::savePCDFileBinary (ss.str(), *pc_right);
        }
        if (pc_rem->size() > 0)
        {
            ss.str("");
            ss << "segments/" << filename << "_rem.pcd";
            std::cout << ss.str() << std::endl;
            pcl::io::savePCDFileBinary(ss.str(), *pc_rem);
        }
    }

    float get_nn_sq_dist(PointType *p)
    {
        if (kdtree.nearestKSearch(*p, 1, idx_search, squared_dist) > 0)
            return squared_dist[0];
        return dist_diag * dist_diag;
    }

    float get_dist_diag()
    {
        return dist_diag / 2;
    }

    float get_dist_apr()
    {
        return dist_apr;
    }
    
    PCDT::Ptr get_cloud()
    {
        return cloud;
    }
    
    PCDT::Ptr get_oct_voxel_cloud()
    {
        return oct_voxel_cloud;
    }
    
    void get_oct_voxel_cloud(const int depth, std::vector<PCDT::Ptr>* z_oct_voxels)
    {
        float z_height = dist_half_edge / std::pow(2,depth - 1);
        if (z_height < p_dist_avg * 2)
            z_height = p_dist_avg * 2;
        int z_bins = (pmax.z - pmin.z) / z_height + 1;
        for (int i = 0; i < z_bins; i++)
        {
            PCDT::Ptr pc(new PCDT);
            z_oct_voxels->push_back(pc);
        }
        
        // filtering for search pcd by octree
        OctreeT octree(dist_half_edge / std::pow(2,depth - 1));
        octree.setInputCloud(cloud);
        octree.defineBoundingBox();
        octree.addPointsFromInputCloud();
        
        Eigen::Vector3f min_pt, max_pt;
        for (auto it = octree.leaf_begin(); it != octree.leaf_end(); it++)
        {
            // get voxel center
            octree.getVoxelBounds(it, min_pt, max_pt);
            PointType p;
            p.x = (min_pt[0]+max_pt[0])/2;
            p.y = (min_pt[1]+max_pt[1])/2;
            p.z = (min_pt[2]+max_pt[2])/2;
            int target_bin = (int)floor((std::min(std::max(p.z, pmin.z), pmax.z) - pmin.z) / z_height);
            // snap to real cloud
            if (kdtree.nearestKSearch(p, 1, idx_search, squared_dist) > 0)
            {
                p.x = (*cloud)[idx_search[0]].x;
                p.y = (*cloud)[idx_search[0]].y;
                p.z = (*cloud)[idx_search[0]].z;
            }
            // form a XYZ cloud
            (*z_oct_voxels)[target_bin]->push_back(p);
        }
        octree.deleteTree();
    }
    
    std::vector<PCDT::Ptr>* get_z_slices()
    {
        return &z_slices;
    }
};

#endif
