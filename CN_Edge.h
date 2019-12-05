#ifndef CN_Edge
#define CN_Edge

#include "CN_Node.h"
#include <vector>
#include <iostream>
#include <unordered_map>
#include <math.h>

using namespace std;

extern void helpRouteOnAtomicEdge(double current, double interest_rate, AtomicEdge* a, Graph* g, int transSeqNum);

// atomic edge, 0 flow, capacity
class AtomicEdge{
public:
	// original Edge Info
	Edge* originEdge;
	int singleCreditIndex;

	bool isDebt;
	bool isLeveraged;
	int atomicEdgeId;
	double capacity;
	double interest_rate;
    int degree;
    bool isAttached;
    double default_rate;
    double collateral_value;

	Node* nodeFrom;
	Node* nodeTo;

	// widget stuff
	WidgetNode* fromWidget;
	WidgetNode* toWidget;


	double AtomicEdge::getSendValue(double IR){
		if (this->isDebt){
			return (1-default_rate) * interest_rate - (1-this->CR()) * default_rate;
		}
		else{
			return collateral_value + interest_rate;
		}
	}

	double AtomicEdge::getReceiveValue(double IR){
		if (this->isDebt){
			return collateral_value + interest_rate;
		}
		else{
			return (1-default_rate) * interest_rate - (1-this->CR()) * default_rate;
		}
	}

	double AtomicEdge::getIR(double value){
		double receiving_ir = (value +  (1-this->CR()) * default_rate) / (1-default_rate);
		// - this->getReceiveValue();
		double sending_ir = value - collateral_value;
		double rr = 0.0;
		if (sending_ir >= receiving_ir){
			return receiving_ir;
		}
		if (receiving_ir > sending_ir){
			return (receiving_ir + sending_ir) / 2;
		}
		return -100;
		// - this->getSendValue();
	}

	void AtomicEdge::updateCollateralValue(double creturns, double areturns){
		double credit_ratio = this->nodeTo->credit_target / (this->nodeTo->credit_target + this->nodeTo->asset_target);
		double return_rate = credit_ratio * creturns + (1-credit_ratio)*areturns;
		if (not this->isDebt){
			double max_extrapolate = this->nodeTo->getWealth(1.0) / this->CR();
			double investment_level = this->nodeTo->credit_target + this->nodeTo->asset_target;
			double max_payment = this->nodeTo->maxCredit(1.0) + this->nodeTo->getScrip();
			double result = 0;
			// double target = min(investment_level,max_payment);
			if (max_extrapolate < investment_level and investment_level < max_payment){
				result = return_rate * (max_extrapolate - investment_level) / max_extrapolate;
			}
			if (investment_level > max_payment){
				result = return_rate * (max_extrapolate - investment_level) / max_extrapolate;
			}			

			// if (max_extrapolate < investment_level and investment_level < max_payment){
				
			// 	result = this->nodeTo->assetReturn * (max_extrapolate - investment_level) / max_extrapolate;
			// 	//negative send
			// }
			// if (max_extrapolate > investment_level and investment_level > max_payment){
			// 	result = this->nodeTo->assetReturn * (max_extrapolate - investment_level) / max_extrapolate;
			// 	//positive send
			// }
			// else {
			// 	result = 0.0;
			// }
			this->collateral_value = result;
		}
		else{
			double hc = 1 + this->CR();
			this->collateral_value = return_rate*(this->nodeTo->maxCredit(hc) - this->nodeTo->maxCredit(1.0));
		}
		// debt, positive receiving
		

	double AtomicEdge::CR(){
		return this->originEdge->singleCreditEdges[singleCreditIndex]->collateralRate;
	}

	AtomicEdge(const AtomicEdge& a, Edge* e, 
			int single, unordered_map<int, AtomicEdge*>& atomicMap, Node* nodeFromT, Node* nodeToT){
		this->originEdge = e;
		this->singleCreditIndex = single;
		this->isDebt = a.isDebt;
		this->atomicEdgeId = a.atomicEdgeId;
		this->capacity = a.capacity;
		this->interest_rate = a.interest_rate;
		atomicMap[this->atomicEdgeId] = this;
		nodeFrom = nodeFromT;
		nodeTo = nodeToT;
		if (this->isDebt){
			nodeFromT->atomicEdge_in[this->atomicEdgeId] = this;
			nodeToT->atomicEdge_out[this->atomicEdgeId] = this;
		} else {
			nodeFromT->atomicEdge_out[this->atomicEdgeId] = this;
			nodeToT->atomicEdge_in[this->atomicEdgeId] = this;
		}
		this->isLeveraged = a.isLeveraged;
		this->isAttached = a.isAttached;
	}

	AtomicEdge(bool d, int id, double capT, double ir, Edge* e, 
		int single, unordered_map<int, AtomicEdge*>& atomicMap, Node* nodeFromT, Node* nodeToT)
		: isDebt(d), atomicEdgeId(id), capacity(capT)
		, interest_rate(ir), originEdge(e), singleCreditIndex(single) {

			atomicMap[this->atomicEdgeId] = this;
			nodeFrom = nodeFromT;
			nodeTo = nodeToT;
			if (this->isDebt){
				nodeFromT->atomicEdge_in[this->atomicEdgeId] = this;
				nodeToT->atomicEdge_out[this->atomicEdgeId] = this;
			} else {
				nodeFromT->atomicEdge_out[this->atomicEdgeId] = this;
				nodeToT->atomicEdge_in[this->atomicEdgeId] = this;
			}
		this->isLeveraged = false;
		this->isAttached = true;		
	}

	void unattach(){
		this->isAttached = false;
	}

	void route(double current, double interest_rate, Graph* g, int transSeqNum){
		if (current == 0.0){
			return;
		}
		// this->print();
		helpRouteOnAtomicEdge(current, interest_rate, this, g, transSeqNum);
	}


	void zero(){
		capacity = 0;
	}

	void print() {
		cout << "    Atomic Edge, id: " << atomicEdgeId << " " 
		<< "capacity: " << capacity << " ir: " << interest_rate 
		<< " isDebt: " << isDebt
		<< " isLeveraged: " << isLeveraged
		<< " single index " << singleCreditIndex 
		<< " from " << nodeFrom->nodeId << " to " << nodeTo->nodeId
		<< endl;
	}

};


class SingleCreditEdge{
public:
	double precision = 100000000;

	double credit_max;
	double credit_interest_rate;
	double collateralRate;
	int maturity;
	bool active;

	AtomicEdge* credit_remain;
	unordered_map<double, AtomicEdge*> debt_current;

	SingleCreditEdge(const SingleCreditEdge& s, Edge* e, 
		int single, unordered_map<int, AtomicEdge*>& atomicMap, Node* nodeFromT, Node* nodeToT) {
		
		this->credit_max = s.credit_max;
		this->credit_interest_rate = s.credit_interest_rate;
		this->collateralRate = s.collateralRate;
		this->active = s.active;
		// this->maturity = s.maturity;
		this->credit_remain = new AtomicEdge(*(s.credit_remain), e, single, atomicMap, nodeFromT, nodeToT);
		for (auto it : s.debt_current){
			this->debt_current[it.first] = 
				new AtomicEdge(*(it.second), e, single, atomicMap, nodeFromT, nodeToT);
		}

	}

	SingleCreditEdge(double c_max, double ir, int& globalId, Edge* e, 
		int single, unordered_map<int, AtomicEdge*>& atomicMap, Node* nodeFromT, Node* nodeToT,double cr,bool active)
		: credit_max(c_max), credit_interest_rate(ir),collateralRate(cr),active(active){
// update atomic edge initiator to include maturity, even for credit edges. may not need maturity for singlecreditedges
			credit_remain = new AtomicEdge(false, 
				globalId++, c_max, ir, e, single, atomicMap, nodeFromT, nodeToT);

	}

 // 	~SingleCreditEdge(){
	// 	delete credit_remain;
	// 	for (auto it : debt_current){
	// 		delete it.second;
	// 	}
	// } 


	void deactivate(){
		this->active = false;
	}

	void activate(){
		this->active = true;
	}	

	void zero(){
		credit_max = 0;
		credit_remain->zero();
	}

	void setCredit(double amt){
		double totalDebt = credit_max - credit_remain->capacity;
		// if(totalDebt<0){
		// 	cout<<"negative debt"<<endl;
		// }
		credit_max = amt;
		credit_remain->capacity = max((amt - totalDebt),0.0);
		
		// cout<<"credit set to "<<amt - totalDebt<<endl;

	}

	double getUtil(){
		double offered = credit_max;
		double used = 0.0;
		for (auto &d : debt_current){
			used += d.second->capacity;
		}
		return max(0.0,used / (offered + 0.00000001));
	}

	bool routeCredit(double current, double interest_rate, int& globalId, Edge* e, 
		int single, unordered_map<int, AtomicEdge*>& atomicMap, Node* nodeFromT, Node* nodeToT){
		if (not this->active){
			cout<<"tried to route on inactive credit line"<<endl;
			return false;
		}
		// cout<<"route credit: "<<current<<"at IR "<<interest_rate<<endl;
		if (round(credit_remain->capacity*precision)/precision < round(precision*current)/precision){
			cout<<"invalid credit capacity"<<endl;
			cout<<"capacity is "<< credit_remain->capacity << " amount is "<< current<<endl;

			return false;
		}

		if (debt_current.find(interest_rate) != debt_current.end()){

			debt_current[interest_rate]->capacity += current;

		} else {
			
			debt_current[interest_rate] = new AtomicEdge(1, 
				globalId++, current, interest_rate, e, single, atomicMap, nodeFromT, nodeToT);

		}

		credit_remain->capacity -= current;
		return true;
	}

	bool routeDebt(double current, double interest_rate) {
		// cout<<"route debt: "<<current<<"at IR "<<interest_rate<<endl;

		if (debt_current.find(interest_rate) == debt_current.end()) {
			cout<<"invalid IR"<<endl;
			return false;
		}
		if (round(precision*debt_current.find(interest_rate)->second->capacity)/precision < round(precision*current)/precision){
			cout<<"invalid debt capacity"<<endl;
			cout<<"capacity is "<< debt_current.find(interest_rate)->second->capacity << " amount is "<< current<<endl;
			return false;
		}

		debt_current.find(interest_rate)->second->capacity -= current;
		if (this->active){			
			credit_remain->capacity = min(current + credit_remain->capacity,credit_max);
		}
		return true;
	}

	void print(){
		cout << "  single credit edge, capacity: " << credit_max 
		<< " ir: " << credit_interest_rate << " cr: "<<collateralRate<<endl;
		credit_remain->print();
		for (auto it : debt_current){
			it.second->print();
		}
	}
};



class Edge{

public:

	vector<int> edgeUsage;
	
	Node* nodeFrom;
	Node* nodeTo;
	vector <SingleCreditEdge*> singleCreditEdges;

	Edge(Node* nodeFromT, Node* nodeToT);
	Edge(const Edge& e, Node* nodeFromT, Node* nodeToT, unordered_map<int, AtomicEdge*>& atomicMap);
	// ~Edge();

	SingleCreditEdge* addSingleCreditEdge(double interest_rate, double cap, 
		int& atomicGlobalId, unordered_map<int, AtomicEdge*>& atomicMap, double collateralRate);

	void routeAtomicEdge(AtomicEdge* a, double flow, double interest_rate, 
		int& atomicGlobalId, unordered_map<int, AtomicEdge*>& atomicMap, int transSeqNum);

	void print();

};


#endif