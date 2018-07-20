#include <iostream>
#include "odas/problems.hpp"
#include "odas/cmaes.hpp"
#include "odas/nsga2.hpp"
#include "odas/popot.hpp"
#include "odas/nlopt.hpp"
#include "odas/voting.hpp"
#include "odas/direct_mlsl.hpp"

int main()
{
    const int seed = 2018;
    Problems probs;
    /** 
    * 9 problems are included for testing purpose :
    * 3 heritage buildings
    * - data/HK_2424275_HKU_Main_Building.pcd"
    * - data/Dublin_25631284_Dublin_City_Hall.pcd"
    * - data/HK_181773526_HKU_HHY_Building.pcd"
    * 3 modern buildings
    * - data/Dublin_177562798_One_Georges_Quay_Plaza.pcd"
    * - data/Dublin_4560539_47_51_OConnell_Street_Upper.pcd"
    * - data/HK_26305856_Fruits_Wholesale_Market.pcd");
    * 3 infrastructures
    * - data/Dublin_46160694_Samuel_Beckett_Bridge.pcd"
    * - data/Dublin_4934685_Sean_OCasey_Bridge.pcd"
    * - data/HK_Two_Piers.pcd
    **/
    for (std::size_t pid = 0; pid < probs.get_problem_size(); pid++)
    {
        std::cout << "problem: " << pid << " " << probs.get_filename(pid) 
                  << std::endl << "============================================"
                  << std::endl;
        {   // GN_DIRECT
            Solvernlopt solver;
            solver.load_prob(probs.get_filename(pid));
            solver.solve(0, 1000);
            solver.validate_best_sol();
        }
        
        {   // compound
            SolverDirectMlsl solver;
            solver.set_seed(seed);
            solver.load_prob(probs.get_filename(pid));
            solver.solve(1000);
            solver.validate_best_sol();
        }
    
        {   // default cmaes 
            Solvercmaes solvercmaes;
            solvercmaes.set_seed(seed);
            solvercmaes.load_prob(probs.get_filename(pid));
            solvercmaes.solve(0, 1000);
            solvercmaes.validate_best_sol();
        }
        {   // sep a IPOP cmaes 
            Solvercmaes solvercmaes;
            solvercmaes.set_seed(seed);
            solvercmaes.load_prob(probs.get_filename(pid));
            solvercmaes.solve(10, 1000);
            solvercmaes.validate_best_sol();
        }
        {   // GN_MLSL_LDS
            Solvernlopt solver;
            solver.set_seed(seed);
            solver.load_prob(probs.get_filename(pid));
            solver.solve(2, 1000);
            solver.validate_best_sol();
        }
        {   // LN_NEWUOA_BOUND
            Solvernlopt solver;
            solver.set_seed(seed);
            solver.load_prob(probs.get_filename(pid));
            solver.solve(4, 1000);
            solver.validate_best_sol();
        }
        {   // Artificial Bee Colony 
            Solverpopot solverabc;
            solverabc.set_seed(seed);
            solverabc.load_prob(probs.get_filename(pid));
            solverabc.solve(2, 1000);
            solverabc.validate_best_sol();
        }
        {   // standard PSO 2007
            Solverpopot solverpso;
            solverpso.set_seed(seed);
            solverpso.load_prob(probs.get_filename(pid));
            solverpso.solve(4, 1000);
            solverpso.validate_best_sol();
        }
        {   // NSGA-II
            Solvernsga2 solverga;
            solverga.set_seed(seed);
            solverga.load_prob(probs.get_filename(pid));
            solverga.solve(1000);
            solverga.validate_best_sol();
        }
        /*
        std::cout << pid << ": ";
        {   // Voting
            Solvervoting solver;
            solver.load_prob(probs.get_filename(pid));
            solver.solve(10);
            solver.validate_best_sol();
        } 
        */
        std::cout << std::endl << std::endl;
    }
    return 0;
}


