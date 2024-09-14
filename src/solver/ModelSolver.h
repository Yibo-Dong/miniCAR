#ifndef MODEL_SOLVER_H
#define MODEL_SOLVER_H

#include "SATInterface.h"
#include "definition.h"
namespace car
{

    /**
     * Categories of flags:
     * 1. O flags used in main solver.      ==> MFlag
     * 2. O flags used in propagation solver.   ==> PFlag
     * 3. propagation flags used in propagation test(to add temporary clauses) ==> PTFlag
     */
    enum SolverFlagEnum
    {
        SolverFlag_Main = 0,
        SolverFlag_Prop = 1,
        SolverFlag_PropTemp = 2,
    };

    class FlagManager
    {
        static int max_flag;

    public:
        static void setStartFlag(int idx) { max_flag = idx; }

        // 1) O flags used in main solver. Used to (de)activate particular frames of O frame.
        static std::vector<int> MFlags;
        // 2) Temporary flags used in inductive checking(propagation).
        static std::vector<int> PTFlags;
        // 3) Propagation flags
        static std::vector<int> PFlags;

    public:
        static int LevelOf(const int flag, SolverFlagEnum ftype);
        static int FlagOf(const int level, SolverFlagEnum ftype);
        static int getNewPTFlag();

    private:
        static int MFlagOf(const int frame_level);
        static int PFlagOf(const int frame_level);
        static int MLevelOf(const int Mflag);
        static int PLevelOf(const int Mflag);
    };

    class ModelSolver
    {
    public:
        // SAT solver
        SATSolver *SATslv;
        // Problem instance, containing TransitionFn and size info.
        const Problem *m;

    private:
        void loadTrans();
        void loadTransSimp();

    public:
        ModelSolver(const Problem *model, bool simp, const SatSolverEnum& sat_option) : m(model)
        {
            SATslv = SATSolverFactory::create_solver(sat_option, model->max_id());
            if (simp)
                loadTransSimp();
            else
                loadTrans();
        }
        ~ModelSolver() { delete SATslv; }

    public:
        /**
         * MFlag / PFlag:
         *      get flag of the corresponding level, together with the cu;
         * PTFlag:
         *      get the negation of latest PTflag, together with the cu;
         * => negated, add to the solver.
         */
        void addNotCubeToLevel(const Cube &cu, int level, SolverFlagEnum ftype, bool primed);

        bool zeroStepReachable(const Cube &from, int target);
        bool oneStepReachable(const Cube &from, int level, bool forwardT);
        void activateLevel(int level, SolverFlagEnum ftype);

        State *getState(bool shrink_to_previous);
        void shrinkModelToPrevious(Assignment &);

        Cube getUCofLatch(bool shrink_to_previous);
    };

}
#endif