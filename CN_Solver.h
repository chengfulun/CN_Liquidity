#ifndef CN_Solver
#define CN_Solver

// #include "CN_WidgetGraph.h"
#include "CN_Constants.h"
#include "CN_CplexConverter.h"
#include <string>
#include <ilcplex/ilocplex.h>
  ILOSTLBEGIN

class LpSolver{
public:
	bool solveLpProblem(CplexConverter& cplexConverter, string mode, string purpose, double cRate, double dRate, double haircut);
	
	void populatebyrow(CplexConverter& cplexConverter, 
		IloModel model, IloNumVarArray x, IloRangeArray c, string mode, string purpose, double cRate, double dRate, double haircut);

	void addObjective(string mode, CplexConverter& cplexConverter, IloModel model, 
	IloNumVarArray x, IloRangeArray c);
};


#endif
