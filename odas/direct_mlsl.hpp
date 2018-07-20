#ifndef __ODAS_DIRECT_MLSL
#define __ODAS_DIRECT_MLSL

#include "types.hpp"
#include "nlopt.hpp"

class SolverDirectMlsl: public Solver
{
    
public:
    SolverDirectMlsl() : Solver("DIRECT_MLSL")
    {
    }
    
    double solve(const int max_eval = 500)
    {
        best_value = 1.0E+8;
        timecost = 0.0;
        //  0 (DIRECT), 2 (GN_MLSL_LDS)
        {
            Solvernlopt solver_direct;
            solver_direct.set_seed(seed);
            solver_direct.load_prob(pcd);
            solver_direct.solve(0, max_eval / 2);
            
            timecost += solver_direct.timecost;
            best_value = solver_direct.best_value;
            best_solution[0] = solver_direct.best_solution[0];
            best_solution[1] = solver_direct.best_solution[1];
        }
        {
            Solvernlopt solver_mlsl;
            solver_mlsl.set_seed(seed);
            solver_mlsl.load_prob(pcd);
            solver_mlsl.solve(2, max_eval / 2);
            
            timecost += solver_mlsl.timecost;
            if (solver_mlsl.best_value < best_value)
            {
                best_value = solver_mlsl.best_value;
                best_solution[0] = solver_mlsl.best_solution[0];
                best_solution[1] = solver_mlsl.best_solution[1];
            }
        }
        
#ifdef TEST_VERBOSE
        std::cerr << "[End] " << name << " : f="<< best_value << ", in " << timecost << "ms." << std::endl;
#endif
        return best_value;
    }
};

#endif
