#ifndef __ODAS_CMAES
#define __ODAS_CMAES

#include "types.hpp"
#include "solver.hpp"

#include "libcmaes/cmaes.h"

class Solvercmaes : public Solver
{
private:
    static const int my_alg_num = 15;
    
public:
    Solvercmaes() : Solver("libcmaes")
    {
    }
    
    static int get_alg_num()
    {
        return my_alg_num;
    }
    
    double solve(const int alg = 10, const int max_eval = 500)
    {
        // good alg: 1 (IPOP_CMAES), 10 (sepaIPOP_CMAES) iter 1000
        const float rho_max = pcd->get_dist_diag();
        max_evaluations = max_eval;
        alg_id = alg;
        
        using namespace libcmaes;
        // Lambda function of evaluation of a solution at client software
        FitFunc fitfunc = [&](const double *x, const int N)
        {
            double ret = evaluate((float)(x[0] * x[0] * rho_max), (float)(x[1] * 2 * M_PI));
            return ret;
        };

        // set up libcmaes
        std::vector<double> x0(dimension, 0.5);
        x0[0] = 0.0;
        x0[1] = 0.5;
        double lbounds[dimension],ubounds[dimension];
        for (int i=0; i<dimension; i++)
        {
            lbounds[i] = 0.0;
            ubounds[i] = 1.0;
        }
        // genotype / phenotype transform associated to bounds.
        GenoPheno<pwqBoundStrategy> gp(lbounds, ubounds, dimension);
        CMAParameters<GenoPheno<pwqBoundStrategy> >
                cmaparams(x0, 0.5, -1, 0, gp); // -1 for automatically decided lambda

        cmaparams.set_seed(seed);
        cmaparams.set_max_fevals(max_eval);
        //cmaparams.set_noisy();
        //cmaparams.set_sep();
        //cmaparams.set_mt_feval(true);
        //cmaparams.set_stopping_criteria(TOLX, false);
        //cmaparams.set_stopping_criteria(STAGNATION,false);
        //cmaparams.set_ftarget(0.0);
        cmaparams.set_elitism(1); // elistism, 0: deactivated, 
                                    // 1: reinjects best seen candidate, 
                                    // 2: initial elitism, reinjects x0, 
                                    // 3: on restart scheme, useful when 
                                    //      optimizer appears to converge to a 
                                    //      value that is higher than the best  
                                    //      value reported along the way
        //if (alg >= 6)
        //    cmaparams.set_tpa(2);   // two-point adapation for step-size update
                                    // 0: no, 1: auto, 2: yes
        cmaparams.set_algo(alg);
        // Available algorithm : CMAES_DEFAULT = 0, IPOP_CMAES, BIPOP_CMAES, aCMAES,
        //                      aIPOP_CMAES = 4, aBIPOP_CMAES, sepCMAES = 6, sepIPOP_CMAES,
        //                      sepBIPOP_CMAES = 8, sepaCMAES, sepaIPOP_CMAES,
        //                      sepaBIPOP_CMAES = 11, VD_CMAES, VD_IPOP_CMAES,
        //                      VD_BIPOP_CMAES = 14

        gettimeofday(&solver_start, NULL);

        // -------------- start solving -------
        cmaes<GenoPheno<pwqBoundStrategy> >(fitfunc, cmaparams);
        // -------------- end of solving

        gettimeofday(&solver_end, NULL);
        long sec  = solver_end.tv_sec  - solver_start.tv_sec;
        long usec = solver_end.tv_usec - solver_start.tv_usec;
        timecost = ((sec) * 1000 + usec/1000.0) + 0.5;
#ifdef TEST_VERBOSE
        std::cerr << "[End] " << name << " (alg "<< alg << "): f="<< best_value << ", in " << timecost << "ms." << std::endl;
#endif
        // checking for unspent budget
        for (int repeat = 1; evaluations < max_evaluations; repeat++)
        {
            Solvercmaes twin;
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
