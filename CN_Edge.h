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
		// nodeFrom is the creditor that holds debt
		if (this->isDebt){
			nodeFromT->atomicEdge_in[this->atomicEdgeId] = this;
			nodeToT->atomicEdge_out[this->atomicEdgeId] = this;
		} 
		// nodeFrom is the creditor that can be paid using debt
		else {
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
		<< "capacity: " << capacity << " ir: " << interest_rate << " cr: " << CR()
		<< " isDebt: " << isDebt
		<< " isLeveraged: " << isLeveraged
		<< " single index " << singleCreditIndex 
		<< " from " << nodeFrom->nodeId << " to " << nodeTo->nodeId
		<< endl;
	}


	double CR();

	double getSendValue(double IR){
		double sendVal;
		if(this->nodeTo->isMarket or this->nodeFrom->isMarket){
			return 0.0;		
		}
		if (this->default_rate > 10){
			this->default_rate = 0.0;
		}
		if (this->collateral_value < -10){
			this->collateral_value = 0.0;
		}

		if (this->collateral_value > 10){
			this->collateral_value = 0.0;
		}

		if (this->isDebt){
			sendVal = -(1-this->default_rate) * IR + (1-this->CR()) * this->default_rate;
			if (this->default_rate > 10000){
				cout<<"dr error"<<endl;
			}
		}
		else{
			sendVal = this->collateral_value - IR;
		}
		if(fabs(sendVal)>10000 ){
			cout<<"sendval too high value"<<endl;
			cout<<"capacity "<<this->capacity<<" nodeFrom: "<<this->nodeFrom->nodeId<<" defaulted "<<this->nodeFrom->defaulted<<" nodeTo: "<<this->nodeTo->nodeId<<" defaulted "<<this->nodeTo->defaulted<<" val: "<<sendVal<<" isDebt: "<<this->isDebt<<" cval: "<<this->collateral_value<<" IR: "<<IR<<" dr: "<<this->default_rate<<endl;
		}

		return max(-0.07,min(0.0,sendVal));
	}

	double getReceiveValue(double IR){
		if(this->nodeTo->isMarket or this->nodeFrom->isMarket){
			return 1000000000000000;
		}
		if (this->default_rate > 10){
			this->default_rate = 0.0;
		}
		if (this->collateral_value < -10){
			this->collateral_value = 0.0;
		}

		if (this->collateral_value > 10){
			this->collateral_value = 0.0;
		}

		double rVal;
		if (this->isDebt){
			rVal = this->collateral_value + IR;
		}
		else{
			rVal = (1-this->default_rate) * IR - (1-this->CR()) * this->default_rate;
		}
		if(fabs(rVal)>10000 ){
			cout<<"receivingVal too high value"<<endl;
			cout<<"capacity "<<this->capacity<<" nodeFrom: "<<this->nodeFrom->nodeId<<" defaulted "<<this->nodeFrom->defaulted<<" nodeTo: "<<this->nodeTo->nodeId<<" defaulted "<<this->nodeTo->defaulted<<" val: "<<rVal<<" isDebt: "<<this->isDebt<<" cval: "<<this->collateral_value<<" IR: "<<IR<<" dr: "<<this->default_rate<<endl;
		}
		return min(0.07,max(0.0,rVal));	
	}

	void updateDR(double DR){
		if(this->nodeTo->isMarket or this->nodeTo->defaulted){
			this->default_rate = 0.0;
		}
		else{
			this->default_rate = DR;
		}
		if(this->default_rate > 1.0){
			cout<<"default rate error"<<endl;
		}
	}

	double getIR(double value){
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

	void updateCollateralValue(double creturns, double areturns){


		double credit_ratio = this->nodeTo->credit_target / (this->nodeTo->credit_target + this->nodeTo->asset_target);
		// cout<<"nodenum "<<this->nodeTo->nodeId<<endl;
		// cout<<"credit_ratio: "<<credit_ratio<<endl;
		// cout<<"credit_target: "<<this->nodeTo->credit_target<<endl;
		// cout<<"asset_target: "<<this->nodeTo->asset_target<<endl;		
		// cout<<"defaulted? "<<this->nodeTo->defaulted<<endl;
		// cout<<"wealth? "<<this->nodeTo->getWealth(1.0)<<endl;
		// cout<<"creturns: "<<creturns<<endl;
		// cout<<"areturns: "<<areturns<<endl;
		double return_rate = credit_ratio * min(0.5,creturns) + (1-credit_ratio)*areturns;
		if(fabs(return_rate > 10)){
			cout<<"error: rr too high! "<<"return rate: "<<return_rate<<" ar: "<<areturns<<" cr: "<<creturns<<endl;
		}

		if(fabs(return_rate > 10)){
			return_rate = areturns;
		}		
		if (not this->isDebt){
			//sum of credit and cash, if all collateral rates were set at the current edge's collateral rate
			// so, sum all credit capacities plus cash if CR = 0
			double max_extrapolate = this->nodeTo->maxExtrapolate(this->CR());
			// getWealth(1.0) / this->CR();
			double investment_level = this->nodeTo->credit_target + this->nodeTo->asset_target;
			double max_payment = this->nodeTo->maxCredit() + this->nodeTo->getScrip();
			double result = 0;
			// double target = min(investment_level,max_payment);
			
			// when taking the collateral requirement hurts an investment that is feasible
			if (max_extrapolate < investment_level and investment_level < max_payment){
				result = return_rate * (max_extrapolate - investment_level) / max_extrapolate;
			}
			// when an investment is infeasible
			if (investment_level > max_payment){
				result = return_rate * (max_extrapolate - max_payment) / max_extrapolate;
			}


			// cout<<"Collateral values check, credit"<<endl;
			// cout<<"rr: "<<return_rate<<" ar: "<<areturns<<" cr: "<<creturns<<" extrapolated investment: "<<max_extrapolate<<" investment level "<<investment_level<<endl;
			
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
			if(result > 10){
				// cout<<"error CV too high!"<<endl;
				result = return_rate;
			}
			if(result < -10){
				result = -return_rate;
			}
			if(result > 0){
				// cout<<"error: why positive sendval?"<<endl;
				// cout<<result<<endl;
			}
			this->collateral_value = result;
		}
		else{
			double cr_bonus = this->CR();
			this->collateral_value = return_rate*(this->nodeTo->maxCredit(cr_bonus,1.0) - this->nodeTo->maxCredit())/(this->nodeTo->maxCredit());
			if(this->nodeTo->maxCredit()==0){
				this->collateral_value = 0;
			}
			if(fabs(this->collateral_value) > 100){
				cout<<"CV too high!"<<endl;
				cout<<"return rate "<<return_rate<<" ar: "<<areturns<<" cr: "<<creturns<<" maxCredit hc "<<this->nodeTo->maxCredit(1.0,cr_bonus)<<" maxcredit "<<this->nodeTo->maxCredit()<<endl;
			}		
		}
		if(this->nodeTo->isMarket){
			this->collateral_value = 0.0;
		}
		if (this->isDebt and this->nodeFrom->defaulted){
			this->collateral_value = 0.0;
		}
		// debt, positive receiving
		if(this->collateral_value > 10){
			this->collateral_value = return_rate;
		}
		if(this->collateral_value < -10){
			this->collateral_value = -return_rate;
		}
	}

	// void printCollateralValues(){

	// }
};




class SingleCreditEdge{
public:
	double precision = 100000000;

	double credit_max;
	double credit_interest_rate;
	double collateralRate;
	int maturity;
	// bool active = true;

	AtomicEdge* credit_remain;
	unordered_map<double, AtomicEdge*> debt_current;

	SingleCreditEdge(const SingleCreditEdge& s, Edge* e, 
		int single, unordered_map<int, AtomicEdge*>& atomicMap, Node* nodeFromT, Node* nodeToT) {
		
		this->credit_max = s.credit_max;
		this->credit_interest_rate = s.credit_interest_rate;
		this->collateralRate = s.collateralRate;
		// this->active = s.active;
		// this->maturity = s.maturity;
		this->credit_remain = new AtomicEdge(*(s.credit_remain), e, single, atomicMap, nodeFromT, nodeToT);
		for (auto it : s.debt_current){
			this->debt_current[it.first] = 
				new AtomicEdge(*(it.second), e, single, atomicMap, nodeFromT, nodeToT);
		}

	}

	SingleCreditEdge(double c_max, double ir, int& globalId, Edge* e, 
		int single, unordered_map<int, AtomicEdge*>& atomicMap, Node* nodeFromT, Node* nodeToT,double cr)
		// ,bool active)
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


	// void deactivate(){
	// 	this->active = false;
	// }

	// void activate(){
	// 	this->active = true;
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
		// if (not this->active){
		// 	cout<<"tried to route on inactive credit line"<<endl;
		// 	return false;
		// }
		if (interest_rate > 0){
			// cout<<"node "<<this->credit_remain->nodeTo->nodeId<<" route credit: "<<current<<" at IR "<<interest_rate<<endl;			
		}
		// this->print();
		if (round(credit_remain->capacity*precision)/precision < round(precision*current)/precision){
			cout<<"error: invalid credit capacity"<<endl;
			cout<<"IR is "<<interest_rate<<"capacity is "<< credit_remain->capacity << " amount is "<< current<<endl;
			this->print();
			return false;
		}

		if (debt_current.find(interest_rate) != debt_current.end()){

			debt_current[interest_rate]->capacity += current;

		} else {
			
			debt_current[interest_rate] = new AtomicEdge(1, 
				globalId++, current, interest_rate, e, single, atomicMap, nodeFromT, nodeToT);

		}
		// cout<<"credit avail: "<<credit_remain->capacity<<endl;
		// credit_remain->capacity -= current;
		// cout<<"credit remaining: "<<credit_remain->capacity<<endl;
		// cout<<"used up: "<<current<<endl;
		return true;
	}

	bool routeDebt(double current, double interest_rate) {
		// cout<<"route debt: "<<current<<"at IR "<<interest_rate<<endl;

		if (debt_current.find(interest_rate) == debt_current.end()) {
			cout<<"error: invalid IR"<<endl;
			return false;
		}
		if (round(precision*debt_current.find(interest_rate)->second->capacity)/precision < round(precision*current)/precision){
			// cout<<"error: invalid debt capacity"<<endl;
			// cout<<"IR is "<<interest_rate<<" capacity is "<< debt_current.find(interest_rate)->second->capacity << " amount is "<< current<<endl;
			// this->print();			
			return false;
		}

		// if (this->active){			
		// }
		// cout<<"debt avail: "<<debt_current.find(interest_rate)->second->capacity<<endl;
		debt_current.find(interest_rate)->second->capacity -= current;		
		credit_remain->capacity = min(current + credit_remain->capacity,credit_max);
		// cout<<"debt remaining: "<<debt_current.find(interest_rate)->second->capacity<<endl;
		// cout<<"used up: "<<current<<endl;

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