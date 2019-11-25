#ifndef CN_Graph
#define CN_Graph

#include "CN_Node.h"
#include "CN_Edge.h"
#include <list>
#include <unordered_map>
#include <fstream>
#include <set>


class Graph{
public:
	int nodeNum;
	int atomicGlobalId;
	int marketId;
	double expected_deposit;
	double expected_asset_return;
	double asset_volatility;
	double CR;
	double deposit_shock;
	double FFR;
	unordered_map<int, Node*> nodes;
	list<Edge*> edges;
	unordered_map<int, AtomicEdge*> atomicEdges;

	vector<double> collateral_values;
	vector<double> interest_values;
	vector<double> default_values;

	unordered_map<string, int> binNum;
	unordered_map<double,double> collateral_valueToMean;
	unordered_map<double,double> interest_valueToMean;
	unordered_map<double,double> default_valueToMean;
	/////////////////////////////////////////////////////////////////////////
	/* Graph basics */
	/////////////////////////////////////////////////////////////////////////
	Graph();
	Graph(int nodeNum, int marketId);
	// Graph(Graph &graphT);
	// Graph& operator=(Graph &graphT);
	// ~Graph();
	void modifyCredit(int s, int t, double credit);
	void deleteEdges();
	void print();
	void printAtomicEdges();
	void printAtomicIouEdges(ofstream& fout);
	void printAvgAtomicIouEdges();
	
	void addMultiEdge(Node* nodeFrom, Node* nodeTo, double credit_ir, double debt_ir, double currDebt, double cap, double cr);
	void updateNodeDegrees();
	void setRoutePreference(vector<string> &v);
	void setThetas(vector<double> &v);
	void setThetas(vector<string> &v);
	void remove(int fid);
	/////////////////////////////////////////////////////////////////////////
	/* Generate Initial Network */
	/////////////////////////////////////////////////////////////////////////
	// void generateTestGraph2();
	// void generateTestGraph3();
	void genTest0Graph(double threshold, int numIR, double cap, double maxCR, double wealth, int marketId);
	void genMarket0Graph(double deposit, double shock, double wealth, double FFR, int marketId, double CR, double EAR, double volatility, bool setT);
	// void genTest1Graph(double threshold, int numIR, int cap);

};


#endif
