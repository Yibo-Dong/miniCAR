#ifndef SAT_INTERFACE_H
#define SAT_INTERFACE_H

#include "definition.h"
// Glucose:
#include "Solver.h"
#include "SolverTypes.h"
// CaDiCaL:
#include "cadical.hpp"
#include <memory>

namespace car
{
    class SATSolver;
    class GlucoseSolver;
    class CaDiCaLSolver;

    class SATSolver
    {
    public:
        virtual ~SATSolver() {}
        virtual void add_clause(const Cube &) = 0;
        virtual void clear_assumption() = 0;
        virtual void push_assumption(const Cube &) = 0;
        virtual bool solve_assumption() = 0;
        virtual Cube get_uc() = 0;
        virtual Cube get_model() = 0;
        virtual void print_assumption(std::ostream &out_stream) = 0;
        virtual void print_clauses(std::ostream &out_stream) = 0;
    };


    class GlucoseSolver : public SATSolver, public Glucose::Solver
    {
    private:
        Glucose::Lit SAT_lit(int id);   // create the Lit used in SAT solver for the id.
        int lit_id(Glucose::Lit) const; // return the id of SAT lit

    public:
        void add_clause(const Cube &) override;      // add a clause
        void clear_assumption() override;            // clear the assumptions
        void push_assumption(const Cube &) override; // push to the assumptions
        bool solve_assumption() override;            // Solve with the assumptions in _assumption.
        Cube get_uc() override;                      // get UC from SAT solver
        Cube get_model() override;                   // get the model from SAT solver
        void print_assumption(std::ostream &out_stream) override;
        void print_clauses(std::ostream &out_stream) override;
    };

    
    class CaDiCaLSolver : public SATSolver, public CaDiCaL::Solver
    {
    private:
        Cube asm_store; // store current assumption, for get_uc()
        int model_size;

    public:
        explicit CaDiCaLSolver(int model_size):model_size(model_size) {}

        void add_clause(const Cube &) override;      // add a clause
        void clear_assumption() override;            // clear the assumptions
        void push_assumption(const Cube &) override; // push to the assumptions
        bool solve_assumption() override;            // Solve with the assumptions in _assumption.
        Cube get_uc() override;                      // get UC from SAT solver
        Cube get_model() override;                   // get the model from SAT solver
        void print_assumption(std::ostream &out_stream) override;
        void print_clauses(std::ostream &out_stream) override;
    };
    
    class SATSolverFactory
    {
        public:
            static SATSolver* create_solver(const SatSolverEnum& sat_option, int model_size){
                switch (sat_option)
                {
                case SAT_CADICAL:
                    return new CaDiCaLSolver(model_size);
                case SAT_GLUCOSE:            
                default:
                    return new GlucoseSolver();
                    break;
                }
            }
            

    };

}

#endif // SAT_INTERFACE_H