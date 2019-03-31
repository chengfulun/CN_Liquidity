#ifndef CN_Node
#define CN_Node

#include <vector>
#include <iostream>
#include <string>
#include <unordered_map>
#include <tuple>

using namespace std;

// forward declaration
class Graph;
class Edge;
class AtomicEdge;
class WidgetNode;

class Node{

private:
	int lag_wealth;
	int lag_deposits;
	int lag_cash;
	int lag_sumAssets;
	bool updateLags();

public:
	int nodeId;
	int transactionNum;

	int srcNum;
	int destNum;
	int successSrc;
	int successDest;
	vector<int> transSeq;

	string routePreference;
	int degree;
	
	bool defaulted;
	bool leveraged;
	// edges
	unordered_map<int, Edge*> edge_out;
	unordered_map<int, Edge*> edge_in;

	// atomic edges
	unordered_map<int, AtomicEdge*> atomicEdge_in;
	unordered_map<int, AtomicEdge*> atomicEdge_out;

	Node(int id);

	// ~Node();

	int getNodeId();
	double getCurrBalance();
	void updateDegree();

	// void remove();
	void makeLeveraged(bool status);
	void makeDefault();
	void printTransSeq();
	void addModification(int transSeqNum);
	double getCollateral(double dRate);
	double getCash();
	double getDebt();
	void makeMarket();
	double getScrip();
	double getWealth(double haircut=0.8);
	double getLeverage();
	bool isMarket;
	// void payIR();

	void getLambda(double E, double sigma_sq, double E_debt, double sigma_sq_debt, double FFR, double limit);
	double theta;
	double deposits;
	double lambda;
	double w_assets;
	double assetSum();
	vector<pair<int, double>> assets;
	double sumAssets();
	void print(); 
	double creditReturn;
	double creditVol;
	void invest();
	double credit_returns_in;
	double credit_vol_in;
	unordered_map<int,double> creditPayIn;
	unordered_map<int,pair<double,double>> creditPayOut;	
	void buyAssets(double amt, int mat);
	double getCredit();
	double folio_volume; 
	int been_defaulted = 0;
	double getOwe();

	double defaultProb;
	double updateDefaultProb();
};

#endif
