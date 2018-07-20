#ifndef __ODAS_NLOPT
#define __ODAS_NLOPT

#include "types.hpp"
#include "solver.hpp"

#include <iomanip>
#include <iostream>
#include <vector>

#include <nlopt.hpp>

static float rho_max_nlopt;
static Solver *solver_nlopt;

static double evaluate_model_nlopt(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
    return solver_nlopt->evaluate((float)(x[0] * x[0] * rho_max_nlopt), (float)(x[1] * 2 * M_PI));
};

class Solvernlopt: public Solver
{
private:
    static const int my_alg_num = 16;
    
    nlopt::algorithm choose_alg(const int alg)
    {
        switch (alg)
        {
        case (0):
            // GOOD Cand
            return nlopt::algorithm::GN_DIRECT;
            break; 
        case (1):
            return nlopt::algorithm::LN_PRAXIS;
            break;
        case (2):
            // GOOD Cand
            return nlopt::algorithm::GN_MLSL_LDS;
            break;
        case (3):
            return nlopt::algorithm::LN_COBYLA;
            break;
        case (4):
            return nlopt::algorithm::LN_NEWUOA_BOUND;
            break;
        case (5):
            return nlopt::algorithm::LN_NELDERMEAD;
            break;
        case (6):
            return nlopt::algorithm::LN_SBPLX;
            break;
        case (7):
            return nlopt::algorithm::LN_AUGLAG;
            break;
        case (8):
            return nlopt::algorithm::LN_BOBYQA;
            break;
        case (9):
            return nlopt::algorithm::GN_DIRECT_L;
            break;
        case (10):
            return nlopt::algorithm::GN_DIRECT_L_RAND;
            break;
        case (11):
            return nlopt::algorithm::GN_DIRECT_NOSCAL;
            break;
        case (12):
            return nlopt::algorithm::GN_DIRECT_L_NOSCAL;
            break;
        case (13):
            return nlopt::algorithm::GN_DIRECT_L_RAND_NOSCAL;
            break;
        case (14):
            return nlopt::algorithm::GN_ORIG_DIRECT;
            break;
        case (15):
            return nlopt::algorithm::GN_ORIG_DIRECT_L;
            break;
        }
        return nlopt::algorithm::LN_NEWUOA_BOUND;
    }
    
public:
    Solvernlopt() : Solver("NLopt")
    {
    }
    
    static int get_alg_num()
    {
        return my_alg_num;
    }
    
    double solve(const int alg = 4, const int max_eval = 500)
    {
        // good alg : 0 (DIRECT), 2 (GN_MLSL_LDS)
        rho_max_nlopt = pcd->get_dist_diag();
        solver_nlopt = (Solver*)this;
        max_evaluations = max_eval;
        alg_id = alg;
    
        nlopt::srand(seed);
        nlopt::opt opt(choose_alg(alg), 2);
        /** possible values: GN_DIRECT = 0, GN_DIRECT_L = 1,  GN_DIRECT_L_RAND,
                            GN_DIRECT_NOSCAL, GN_DIRECT_L_NOSCAL, 
                            GN_DIRECT_L_RAND_NOSCAL, GN_ORIG_DIRECT,
                            GN_ORIG_DIRECT_L, ***GD_STOGO, ***GD_STOGO_RAND,
                            ***LD_LBFGS_NOCEDAL, ***LD_LBFGS, LN_PRAXIS = 12, 
                            ***LD_VAR1, ***LD_VAR2, ***LD_TNEWTON, 
                            ***LD_TNEWTON_RESTART, ***LD_TNEWTON_PRECOND,
                            ***LD_TNEWTON_PRECOND_RESTART, ***GN_CRS2_LM, 
                            GN_MLSL = 20, GD_MLSL = 21, GN_MLSL_LDS = 22,
                            GD_MLSL_LDS, LD_MMA, LN_COBYLA = 25, LN_NEWUOA = 26,
                            LN_NEWUOA_BOUND = 27, LN_NELDERMEAD = 28, LN_SBPLX = 29, LN_AUGLAG = 30,
                            ***LD_AUGLAG, LN_AUGLAG_EQ = 32, ***LD_AUGLAG_EQ, LN_BOBYQA = 34,
                            ***GN_ISRES, ***AUGLAG, ***AUGLAG_EQ, ***G_MLSL, ***G_MLSL_LDS,
                            ***LD_SLSQP, ***LD_CCSAQ, ***GN_ESCH,
            *** = not DFO
        */
        //std::vector<double> lb(dimension), ub(dimension);
        //lb[0] = 0.0; 
        //lb[1] = 0.0;
        //ub[0] = 1.0; 
        //ub[1] = 1.0;
        opt.set_lower_bounds(0.0);
        opt.set_upper_bounds(1.0);
        opt.set_maxeval(max_evaluations);
        opt.set_min_objective(evaluate_model_nlopt, NULL);
        opt.set_xtol_rel(1e-6);
        std::vector<double> x(dimension);
        x[0] = 0.0; 
        x[1] = 0.5;
        double minf;

        gettimeofday(&solver_start, NULL);
        // -------------- start solving -------
        try
        {
            nlopt::result result = opt.optimize(x, minf);
            //std::cout << "found minimum at f(" << x[0] << "," << x[1] << ") = "
            //    << std::setprecision(10) << minf << std::endl;
        }
        catch(std::exception &e) 
        {
            std::cout << "nlopt failed: " << e.what() << std::endl;
        }
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
            Solvernlopt twin;
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
