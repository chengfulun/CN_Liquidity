#ifndef CN_CplexConverter
#define CN_CplexConverter

// #include "CN_WidgetGraph.h"
#include "CN_Constants.h"
#include <unordered_map>
#include "CN_Node.h"
#include "CN_Graph.h"
class CplexConverter{

public:
	int transSeqNum;

	struct Variable{
		int varId;
		int atomicEdgeId;
		double interest_rate;
		double value_receiving;
		double value_paying;
		Variable(int v, int a, double val_in, double val_out, double IR) : varId(v), atomicEdgeId(a), value_receiving(val_in), value_paying(val_out), interest_rate(IR) {}
	};

	vector<Variable> variables;
	unordered_map <int, vector<int>> atomicIdToVarIdDict;
	unordered_map <int, int> varIdToAtomicIdDict;
	vector<double> capacities;
	vector<double> results;
	double request;
	Node* src;
	Node* dest;
	Graph* graph;
	bool verb;
	
	void constructCplex(Graph* g, Node* s, Node* t, double req, int transSeqNumT, bool verb);

	void printInput();

	void printResult();

	double copyBack();

	double getRouted();

	double getIRCost();

	double getDebtCost();

	double valueCluster(int val, int type);

	double getAssetCost();

	double check();

};


#endif