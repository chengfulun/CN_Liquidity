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

	Node* nodeFrom;
	Node* nodeTo;

	// widget stuff
	WidgetNode* fromWidget;
	WidgetNode* toWidget;

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

	AtomicEdge* credit_remain;
	unordered_map<double, AtomicEdge*> debt_current;

	SingleCreditEdge(const SingleCreditEdge& s, Edge* e, 
		int single, unordered_map<int, AtomicEdge*>& atomicMap, Node* nodeFromT, Node* nodeToT) {
		
		this->credit_max = s.credit_max;
		this->credit_interest_rate = s.credit_interest_rate;
		this->collateralRate = s.collateralRate;
		// this->maturity = s.maturity;
		this->credit_remain = new AtomicEdge(*(s.credit_remain), e, single, atomicMap, nodeFromT, nodeToT);
		for (auto it : s.debt_current){
			this->debt_current[it.first] = 
				new AtomicEdge(*(it.second), e, single, atomicMap, nodeFromT, nodeToT);
		}

	}

	SingleCreditEdge(double c_max, double ir, int& globalId, Edge* e, 
		int single, unordered_map<int, AtomicEdge*>& atomicMap, Node* nodeFromT, Node* nodeToT,double cr)
		: credit_max(c_max), credit_interest_rate(ir),collateralRate(cr){
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

	bool routeCredit(double current, double interest_rate, int& globalId, Edge* e, 
		int single, unordered_map<int, AtomicEdge*>& atomicMap, Node* nodeFromT, Node* nodeToT){

		// cout<<"route credit: "<<current<<"at IR "<<interest_rate<<endl;
		if (round(credit_remain->capacity*precision)/precision < round(precision*current)/precision){
			// cout<<"invalid credit capacity"<<endl;
			// cout<<"capacity is "<< credit_remain->capacity << " amount is "<< current<<endl;

			// return false;
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
			// cout<<"invalid debt capacity"<<endl;
			// cout<<"capacity is "<< debt_current.find(interest_rate)->second->capacity << " amount is "<< current<<endl;
			// return false;
		}

		debt_current.find(interest_rate)->second->capacity -= current;
		credit_remain->capacity = min(current + credit_remain->capacity,credit_max);
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