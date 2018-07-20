#ifndef __ODAS_NSGA2
#define __ODAS_NSGA2

#include "types.hpp"
#include "solver.hpp"

#include <nsga2/global.h>
#include <nsga2/NSGA2.h>

static float rho_max_nsga2;
static Solver *solver_nsga2;

void evaluate_model_nsga2(double *xreal, double *xbin, int **gene, double *obj, double *constr) 
{
    obj[0] = solver_nsga2->evaluate((float)(xreal[0] * xreal[0] * rho_max_nsga2), (float)(xreal[1] * 2 * M_PI));
};

//void no_report_nsga2(nsga2::population& pop){};

class Solvernsga2: public Solver
{
public:
    Solvernsga2() : Solver("nsga2")
    {
    }
    
    double solve(const int max_eval = 500)
    {
        rho_max_nsga2 = pcd->get_dist_diag();
        solver_nsga2 = (Solver*)this;
        max_evaluations = max_eval;
    
        nsga2::individual_config conf;
        conf.limits_realvar.push_back(std::make_pair(0.0,1.0));
        conf.limits_realvar.push_back(std::make_pair(0.0,1.0));
        
        nsga2::NSGA2 nsga2;
        nsga2.set_seed(seed);
        nsga2.set_nreal(2);
        nsga2.set_nbin(0);
        nsga2.set_nobj(1);
        nsga2.set_ncon(0); // constraint 
        nsga2.set_popsize(40);
        nsga2.set_ngen(max_evaluations*2/40);
        nsga2.set_pcross_real(0.9);
        nsga2.set_pcross_bin(0.0);
        nsga2.set_pmut_real(0.1);
        nsga2.set_pmut_bin(0.0);
        nsga2.set_eta_c(10);
        nsga2.set_eta_m(10);
        nsga2.set_epsilon_c(1e-14);
        nsga2.set_limits_realvar(conf.limits_realvar);
        nsga2.set_function(evaluate_model_nsga2);
        nsga2.set_crowdobj(false); // crowd over the parameters, not the objective functions
        //nsga2.set_crowdobj(true); // crowd over objective function
        // nsga2.set_custom_report_function(&no_report_nsga2);
        // nsga2.set_nreport(10);
        nsga2.set_backup_filename(""); // no backup
        
        nsga2.initialize();
        
        gettimeofday(&solver_start, NULL);
        // -------------- start solving -------
        nsga2.evolve();
        // -------------- end of solving

        gettimeofday(&solver_end, NULL);
        long sec  = solver_end.tv_sec  - solver_start.tv_sec;
        long usec = solver_end.tv_usec - solver_start.tv_usec;
        timecost = ((sec) * 1000 + usec/1000.0) + 0.5;
#ifdef TEST_VERBOSE
        std::cerr << "[End] " << name << ": f="<< best_value << ", in " << timecost << "ms." << std::endl;
#endif
        // checking for unspent budget
        for (int repeat = 1; evaluations < max_evaluations; repeat++)
        {
            Solvernsga2 twin;
            twin.set_seed(seed + repeat * 2018);
            twin.load_prob(pcd);
            twin.solve(max_evaluations - evaluations);
            
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
