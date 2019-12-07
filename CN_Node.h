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
	void makeLeveraged(bool status, double haircut);
	void makeDefault();
	void printTransSeq();
	void addModification(int transSeqNum);
	double getCollateral(double dRate);
	double getCash();
	double getDebt();
	void makeMarket();
	double getScrip();
	double getWealth(double haircut);
	double credit2Return(double DR);
	void postCashUpdate(double haircut);
	bool isMarket;
	// void payIR();

	vector<double> credit_investments;
	vector<double> credit_returns; // aggregate revenue

	void getLambda(double E, double sigma_sq, double E_debt, double sigma_sq_debt, double FFR, double mlimit, double haircut);
	double theta;
	double deposits;
	double lambda;
	double w_assets;
	double mReserved;
	double reserveT;

	double assetSum();
	// maturity, amount
	vector<pair<int, double>> assets;
	double sumAssets();
	void print(); 
	// double creditReturn;
	double creditVol;
	void invest();
	double credit_returns_in;
	double credit_vol_in;
	unordered_map<int,double> creditPayIn;
	unordered_map<int,pair<double,double>> creditPayOut;	
	void buyAssets(double amt);
	double getCredit();
	double folio_volume; 
	int been_defaulted = 0;
	double getOwe();
	double getDemand(double price, double assetR);
	double expected_creditrevenue();
	double creditReturn(double DR);

	double credit_target; //executed amount
	double asset_target;  //executed amount
	

	double getInterestOwed();

	// collateral value
	double maxCredit(double haircut);

	void updateWealthReturn();
	void updateAssetReturn();
	void updateCreditReturn();
	void updateAllReturns();
	void updateInterestPaid();

	void updateAllHist(double price, double assetR, double creditR,double creditI,double assetI, double IR, int length);
	// double updateWealthHist();
	// double updateAssetHist();
	// double updateCreditHist();

	vector<double> wealth_history;
	vector<double> asset_returns;
	vector<double> asset_investments;
	vector<double> asset_leverage_interest;

	double wealth_return;
	double asset_return;
	double credit_return;
	double investment_return;
	double interest_paid;

	double asset_aggression;
	double credit_aggression;

	double creditRequestTarget = 0.0;


};

#endif