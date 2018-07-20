#ifndef __ODAS_VOTING
#define __ODAS_VOTING

#include "types.hpp"
#include "solver.hpp"


class Solvervoting : public Solver
{
private:
    static const int FULL_PCD = 0;
    static const int OCTREE1 = 1;
    static const int OCTREE2 = 2;
    static const int OCTREE3 = 3;
    static const int OCTREE4 = 4;
    static const int OCTREE5 = 5;
    static const int OCTREE6 = 6;
    static const int OCTREE7 = 7;
    static const int OCTREE8 = 8;
    static const int OCTREE9 = 9;
    static const int Z_SLICES = 10;
    static const int SIFT_PCD = 11;
    
    const int num_alg = 9;
    int rho_resolution = 100;
    int theta_resolution = 1000;
    const float M_HALF_PI = M_PI / 2;
    const float M_ONE_AND_HALF_PI = M_PI * 3 / 2;
    const float M_TWO_PI = M_PI * 2;
public:
    Solvervoting() : Solver("voting")
    {
    }
    
    int get_alg_num()
    {
        return num_alg;
    }
    
    long long int get_key(float x1, float y1, float x2, float y2)
    {
        float mx = (x1 + x2) / 2;
        float my = (y1 + y2) / 2;
        float l_theta = std::fmod(atan2f(y2 - y1, x2 - x1) + M_PI/2, M_PI);;
        float b = my - mx * tan(l_theta);
        float rho = fabsf(b * std::cos(l_theta));
        float theta = l_theta;
        if (b >= 0)
        {
            if (l_theta <= M_HALF_PI)
                theta += M_HALF_PI;
            else
                theta -= M_HALF_PI;
        }
        else
        {
            if (l_theta <= M_HALF_PI)
                theta -= M_HALF_PI;
            else
                theta -= M_ONE_AND_HALF_PI;
        }
        if (theta > M_TWO_PI)
            theta -= M_TWO_PI;
        if (theta < 0)
            theta += M_TWO_PI;
        if (theta > M_TWO_PI)
            theta -= M_TWO_PI;
        if (theta < 0)
            theta += M_TWO_PI;
        long long ret = (long long)std::round(rho * rho_resolution);
        ret *= 10 * theta_resolution;
        ret += (long long)std::round(theta / M_PI * theta_resolution);
        return ret;
    }
    
    float key_to_rho(long long int key)
    {
        int long_rho = (int)(key / (10 * theta_resolution));
        return long_rho * 1.0f / rho_resolution;
    }
    
    float key_to_theta(long long int key)
    {
        int long_theta = (int)(key % (10 * theta_resolution));
        return long_theta * M_PI / theta_resolution;
    }
    
    double solve(const int alg = 1)
    {
        alg_id = alg;
        const float rho_max = pcd->get_dist_diag();
        const float dist_apr = pcd->get_dist_apr();
        
        std::vector<PCDT::Ptr>* parts = nullptr;
        if (alg == Z_SLICES)
        {
            parts = pcd->get_z_slices();
            rho_resolution = 100;
            theta_resolution = 200;
        }
        else
        {
            parts = new std::vector<PCDT::Ptr>();
            switch (alg)
            {
            case (FULL_PCD):
                //PCDT::Ptr pc(pcd->get_cloud());
                //parts->push_back(pc);
                rho_resolution = 1000;
                theta_resolution = 1000;
                break;
            case (OCTREE1):
                pcd->get_oct_voxel_cloud(1, parts);
                rho_resolution = 10;
                theta_resolution = 10;
                break;
            case (OCTREE2):
                pcd->get_oct_voxel_cloud(2, parts);
                rho_resolution = 30;
                theta_resolution = 10;
                break;
            case (OCTREE3):
                pcd->get_oct_voxel_cloud(3, parts);
                rho_resolution = 100;
                theta_resolution = 10;
                break;
            case (OCTREE4):
                pcd->get_oct_voxel_cloud(4, parts);
                rho_resolution = 100;
                theta_resolution = 30;
                break;
            case (OCTREE5):
                pcd->get_oct_voxel_cloud(5, parts);
                rho_resolution = 100;
                theta_resolution = 100;
                break;
            case (OCTREE6):
                pcd->get_oct_voxel_cloud(6, parts);
                rho_resolution = 100;
                theta_resolution = 300;
                break;
            case (OCTREE7):
                pcd->get_oct_voxel_cloud(7, parts);
                rho_resolution = 100;
                theta_resolution = 300;
                break;
            case (OCTREE8):
                pcd->get_oct_voxel_cloud(8, parts);
                rho_resolution = 200;
                theta_resolution = 300;
                break;
            case (OCTREE9):
                pcd->get_oct_voxel_cloud(9, parts);
                rho_resolution = 200;
                theta_resolution = 300;
                break;
            }
        }
        std::unordered_map<long long int, int> votes;
        long long int key, max_key;
        int max_cnt = 0;
        
        gettimeofday(&solver_start, NULL);
        // -------------- start solving -------
        for(auto pc : *parts)
        {
            auto pc_end = pc->end();
            // each slides of pc
            for (auto it = pc->begin(); it != pc_end; it++)
            {
                for (auto iter = it + 1; iter != pc_end; iter++)
                {
                    if (fabs(it->z - iter->z) > dist_apr)
                        continue;
                    key = get_key(it->x, it->y, iter->x, iter->y);
                    auto result = votes.find(key);
                    if (result == votes.end())
                    {
                        // not found
                        votes.insert(std::make_pair(key, 1));
                    }
                    else
                    {
                        // found
                        result->second ++;
                        if (result->second > max_cnt)
                        {
                            max_key = key;
                            max_cnt = result->second;
                        }
                    }
                }   // for each point pair
            }   // for each point
        }   // for each parts
        
        best_solution[0] = key_to_rho(max_key);
        best_solution[1] = key_to_theta(max_key);
        // -------------- end of solving ------
        gettimeofday(&solver_end, NULL);
        long sec  = solver_end.tv_sec  - solver_start.tv_sec;
        long usec = solver_end.tv_usec - solver_start.tv_usec;
        timecost = ((sec) * 1000 + usec/1000.0) + 0.5;
#ifdef TEST_VERBOSE
        std::cerr << "[End] " << name << ": f="<< best_value << ", in " << (((sec) * 1000 + usec/1000.0) + 0.5) << "ms." << std::endl;
#endif
        return best_value;
    }
};

#endif
