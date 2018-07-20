#ifndef __ODAS_SOLVER
#define __ODAS_SOLVER

#include "types.hpp"
#include "simpcd.hpp"

class Solver
{
private:
    static const int alg_num = 1;
    bool protected_pcd;
    
protected:
    const int dimension = 2;
    int evaluations, max_evaluations;
    int seed, alg_id;
    std::string name;
    SymPCD* pcd;
    struct timeval solver_start, solver_end;
    
public:
    double best_value, timecost;
    double valid_apr, valid_rmse;
    float best_solution[2];

    Solver(std::string name_str = "no name")
        : name(name_str)
    {
        best_value = 1.0E+8;
        for (int i = 0; i < dimension; i++)
            best_solution[i] = 0.0f;
        seed = 1;
        alg_id = 0;
        max_evaluations = 500;
        protected_pcd = false;
    }

    ~Solver()
    {
        if (!protected_pcd)
            delete pcd;
    }
    
    static int get_alg_num()
    {
        return alg_num;
    }
    
    void set_seed(int rand_seed)
    {
        seed = rand_seed;
    }

    void load_prob(std::string filename, const int depth = 4, const float epsilon = 0.005)
    {
        pcd = new SymPCD();
        pcd->load(filename, depth, epsilon);
        evaluations = 0;
    }
    
    void load_prob(SymPCD* extern_pcd)
    {
        pcd = extern_pcd;
        protected_pcd = true;
        evaluations = 0;
    }

    double evaluate(const float rho, const float theta)
    {
        if (evaluations >= max_evaluations)
            return 1.0E+8;
        double ret = pcd->fast_approx_fitness(rho, theta);
        if (ret < best_value)
        {
            best_value = ret;
            best_solution[0] = rho;
            best_solution[1] = theta;
            //std::cout << (evaluations+1) << ", ";
            //validate_best_sol();
#ifdef TEST_VERBOSE
            gettimeofday(&solver_end, NULL);
            long sec  = solver_end.tv_sec  - solver_start.tv_sec;
            long usec = solver_end.tv_usec - solver_start.tv_usec;
            std::cerr << "> " << evaluations <<  ", " << (((sec) * 1000 + usec/1000.0) + 0.5) << ", " << rho << ", "<< theta << ", " << ret << std::endl;
#endif
        }
        ++evaluations;
        return ret;
    }

    double validate_best_sol()
    {
        valid_rmse = 100.0f;
        valid_apr = 0.0f;
        if (!std::isnan(best_solution[0]) && !std::isnan(best_solution[1]))
            valid_rmse = pcd->validate(best_solution[0], best_solution[1], &valid_apr);
#ifdef TEST_BENCHMARK
        std::cout << name << "_" << alg_id << ", "<< std::round(timecost)/1000 << "s, f=" << best_value/(pcd->get_dist_diag()*2) << ", apr=" << std::round(valid_apr*10000)/100 << "%, rmse=" << std::round(valid_rmse*1000)/1000 << "m, \n >> Optimal solution X = (" << best_solution[0] << "m, " << best_solution[1]/M_PI << "Pi)" <<  std::endl;
#endif
        return valid_rmse;
    }

    void segment_pcd()
    {
        if (!std::isnan(best_solution[0]) && !std::isnan(best_solution[1]))
        {
            pcd->segment(best_solution[0], best_solution[1]);
        }
    }
};

#endif
