#ifndef __ODAS_POPOT
#define __ODAS_POPOT

#include "types.hpp"
#include "solver.hpp"

#include <jsoncpp/json/json.h>
#include <popot/rng_generators.h>
typedef popot::rng::CRNG RNG_GENERATOR;
#include <popot/popot.h>
typedef popot::Vector<double> TVector;

class Solverpopot: public Solver
{
private:
    static const int my_alg_num = 5;
    
public:
    Solverpopot() : Solver("popot")
    {
    }
    
    static int get_alg_num()
    {
        return my_alg_num;
    }
    
    double solve(const int alg = 4, const int max_eval = 500)
    {
        // good alg cand = 4 (spso2007)
        const float rho_max = pcd->get_dist_diag();
        max_evaluations = max_eval;
        alg_id = alg;
    
        RNG_GENERATOR::rng_srand(seed);
        RNG_GENERATOR::rng_warm_up();

        auto lbound = [] (size_t index) -> double { return 0.0; };
        auto ubound = [] (size_t index) -> double { return 1.0; };
        auto stop =   [&] (double fitness, int epoch) -> bool { return evaluations >= max_evaluations;};
        auto cost_function = [&] (TVector &pos) -> double { 
            return evaluate((float)pos.getValuesPtr()[0] * pos.getValuesPtr()[0] * rho_max, (float)(pos.getValuesPtr()[1] * 2 * M_PI));
        };
        
        size_t swarm_size = 20; 
        size_t colony_size = 20; 
        size_t nbcats = 10;
        size_t nbloups = 50;
        popot::algorithm::Base* algo;
        gettimeofday(&solver_start, NULL);
        switch (alg)
        {
        case (0):
            algo = popot::algorithm::gwo(nbloups, dimension, lbound, ubound, max_evaluations, stop, cost_function);
            break;
        case (1):
            algo = popot::algorithm::cso(nbcats, dimension, lbound, ubound, stop, cost_function);
            break;
        case (2):
            // Artificial Bee Colony 
            algo = popot::algorithm::abc(colony_size, dimension, lbound, ubound, stop, cost_function);
            break;
        case (3):
            algo = popot::algorithm::spso2006(dimension, lbound, ubound, stop, cost_function);
            break;
        case (4):
            algo = popot::algorithm::spso2007(dimension, lbound, ubound, stop, cost_function);
            break;
        default:
            algo = popot::algorithm::abc(colony_size, dimension, lbound, ubound, stop, cost_function);
            break;
        }
        // -------------- start solving -------
        algo->init();
        algo->run();
        // -------------- end of solving ------
        gettimeofday(&solver_end, NULL);
        delete algo;

        long sec  = solver_end.tv_sec  - solver_start.tv_sec;
        long usec = solver_end.tv_usec - solver_start.tv_usec;
        timecost = ((sec) * 1000 + usec/1000.0) + 0.5;
#ifdef TEST_VERBOSE
        std::cerr << "[End] " << name << " (alg "<< alg << "): f="<< best_value << ", in " << timecost << "ms." << std::endl;
#endif
        // checking for unspent budget
        for (int repeat = 1; evaluations < max_evaluations; repeat++)
        {
            Solverpopot twin;
            twin.set_seed(seed + repeat * 2018);
            twin.load_prob(pcd);
            twin.solve(alg, max_evaluations - evaluations);
            
            evaluations += twin.evaluations;
            timecost += twin.timecost;
            if (best_value > twin.best_value)
            {
                best_value = twin.best_value;
                best_solution[0] = twin.best_solution[0];
                best_solution[1] = twin.best_solution[1];
            }
        }
        return best_value;
    }
};

#endif
