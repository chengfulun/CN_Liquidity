// Node.C

#include "CN_Node.h"
#include "CN_Edge.h"
#include <queue>
#include <vector>
#include <fstream>
#include <math.h>
#include <random>
#include <chrono>
#include <algorithm>

using namespace std;

extern void printEdge(Edge*);
extern void printAtomicEdge(AtomicEdge*);

Node::Node(int id){
	this->nodeId = id;
	this->transactionNum = 0;
	this->routePreference = "";
	this->srcNum = 0;
	this->destNum = 0;
	this->successSrc = 0;
	this->successDest = 0;
	// this->defaulted = false;
	this->leveraged = false;
	this->isMarket = false;
	this->creditReturn = 0;
	this->creditVol = 0;
	this->folio_volume = 0.0;

	this->lag_wealth = -1.0;
	this->lag_deposits = -1.0;
	this->lag_cash = -1.0;
	this->lag_sumAssets = -1.0;
}


// Node::~Node(){
// 	for (auto it : edge_in){
// 		delete it.second;
// 	}
// 	for (auto it : edge_out){
// 		delete it.second;
// 	}
// 	for (auto it : atomicEdge_in){
// 		delete it.second;
// 	}
// 	for (auto it : atomicEdge_out){
// 		delete it.second;
// 	}	
// }

void Node::makeMarket(){
	this->isMarket = true;
}

void Node::makeDefault(){
	// for (auto &itIn : atomicEdge_in){
	// 	itIn.second->zero();
	// }	
	// for (auto &itOut : atomicEdge_out){
	// 	itOut.second->zero();
	// }
	this->defaulted = true;
	// this->creditPayOut.clear();
}

/**
 * if status is True: check if each individual edge is over-leveraged
 * if status is False: update so no edge is over-leveraged
 */
void Node::makeLeveraged(bool status){
	if (status){
		this->leveraged = true;
		double ratio = this->getWealth(0.8) / (this->getDebt() + deposits);
		for (auto &it : atomicEdge_in){
			if (not it.second->isDebt){
				int cIx = it.second->singleCreditIndex;
				double cRate = it.second->originEdge->singleCreditEdges[cIx]->collateralRate;
				if(ratio < cRate && not this->isMarket){
					it.second->isLeveraged = true;
				}
			}
		}
	}
	else{
		this->leveraged = false;
		for (auto &it : atomicEdge_in){
			if (not it.second->isDebt){
				it.second->isLeveraged = false;				
			}
		}		
	}
}

/**
 * Calculates the leverage
 * TODO do not hard code haircut=0.8
 */
double Node::getLeverage(){
	return this->getWealth(0.8) / (this->getDebt() + deposits);
}

void Node::print(){
	cout << "Node " << this->nodeId << endl;
	for (auto &it : edge_in){
		printEdge(it.second);
	}
	for (auto &it : edge_out){
		printEdge(it.second);
	}
	cout << "Atomic Edges: " << endl;
	for (auto &it : atomicEdge_in){
		printAtomicEdge(it.second);
	}
	for (auto &it : atomicEdge_out){
		printAtomicEdge(it.second);
	}
	cout << endl;

	return;
}

void Node::updateDegree(){
	this->degree = edge_in.size();
}

int Node::getNodeId(){ return nodeId; }

double Node::getCurrBalance(){
	double temp = 0;
	for (auto& it : this->atomicEdge_in){
		if (it.second->isDebt){
			temp += it.second->capacity * it.second->interest_rate;
		}
	}
	for (auto& it : this->atomicEdge_out){
		if (it.second->isDebt){
			temp -= it.second->capacity * it.second->interest_rate;
		}
	}
	return temp;
}

// extern double getNodeCurrBalance(Node*);

// double Node::getCurrBalance(){
// 	return getNodeCurrBalance(this);
// }

void Node::addModification(int transSeqNum){
	this->transSeq.push_back(transSeqNum);
}

void Node::printTransSeq(){
	for (int i = 0; i < transSeq.size(); ++i){
		cout << transSeq[i] << " ";
	}
	cout << endl;
}

// void Node::remove(){
	// for (auto ei : edge_in){
		// delete ei.second;
		// for(auto sei : ei.second->singleCreditEdges){
			// delete sei;
		// }
	// }
// }

double Node::getCollateral(double dRate){
	double temp = 0;
	double amt = 0;
	for (auto &eIn : edge_in){
		for (auto &ceIn : eIn.second->singleCreditEdges){
			double t = (ceIn->credit_max - ceIn->credit_remain->capacity);
			amt += t;
			temp += t*ceIn->collateralRate;
			// cout<<"debt: "<<t<<"at crate "<<ceIn->collateralRate<<endl;
		}
	}
	temp += deposits * dRate;
	// cout<<this->getNodeId()<<" total debt: "<<amt<<endl;
	// cout<<this->getNodeId()<<" total required collateral: "<<temp<<endl;
	return temp;
}

double Node::getCash(){
	if(defaulted){
		return 0.0;
	}
	double temp = 0.0;
	// cout<<nodeNum<<endl;
	for (auto &eOut : edge_out){
		if (eOut.second->nodeTo->isMarket){
			for(auto &it: eOut.second->singleCreditEdges){
				for(auto &d: it->debt_current)
					temp += d.second->capacity;					
			}
		}
	}
	return temp;
}

double Node::sumAssets(){
	double s = 0.0;
	for (auto &a : assets){
		s += a.second;
	}
	return s;
}

double Node::getOwe(){
	double paySum = 0;
	for (auto &toPay : this->creditPayOut){
		paySum += toPay.second.second;
	}	
	return paySum;
}

double Node::getWealth(double haircut){
	if(defaulted){
		return 0.0;
	}
	// if(v>0){
	// // cout<<"node "<<nodeId<<" IOUs "<<this->getScrip()<< " debt "<<this->getDebt()<<" assets "<<this->sumAssets()<<" deposits "<<this->deposits<<" credit "<<getCredit()<<endl;		
	// }
	return -this->getOwe() + this->getScrip() - this->getDebt() + this->sumAssets() * haircut - this->deposits;
	// double temp = 0;
	// cout<<nodeNum<<endl;
}

double Node::getScrip(){
	if(defaulted){
		return 0.0;
	}
	double temp = 0.0;
	for (auto &dOut : atomicEdge_in){
		if (dOut.second->isDebt){
			temp += dOut.second->capacity;
		}
	}
	return temp;
}

/**
 * Calculate total incoming credit (0 of default)
 */
double Node::getCredit(){
	if(defaulted){
		return 0.0;
	}
	double temp = 0.0;
	// cout<<nodeNum<<endl;
	for (auto &eIn : edge_in){
		for(auto &it: eIn.second->singleCreditEdges){
			temp += it->credit_remain->capacity;
		}
	}
	return temp;
}

double Node::getDebt(){
	double temp = 0;
	for (auto &dOut : atomicEdge_out){
		if (dOut.second->isDebt && not dOut.second->nodeFrom->defaulted){
			temp += dOut.second->capacity;
		}
	}
	return temp;
}

/**
 * Calculates three variables described in section 3.4 of AI3 paper
 * lambda: proportion of capital invested that maximizes utility ,
 * w_assets: the optimal proportion of investment allocated to assets
 * folio_volume: actual maximum portfolio size (limited by leverage)
 */
void Node::getLambda(double E, double sigma_sq, double E_debt, double sigma_sq_debt, double FFR, double mlimit){
	// double E_debt = this->credit_returns_in;
	// double sigma_sq_debt = this->credit_vol_in;
	// if(this->leveraged){
	// 	// cout<<"node "<<nodeId<<" leveraged"<<endl;
	// 	this->lambda = 0.0;
	// 	this->w_assets = 0.0;
	// 	this->folio_volume = 0.0;
	// 	return;
	// }
	double floor = 0.0;
	if (E_debt<0){
		E_debt = 0.0;
		std::default_random_engine generator;
		// generator(std::chrono::system_clock::now().time_since_epoch().count());
		std::uniform_real_distribution<double> distribution(0.0,0.2);		
	    floor = distribution(generator);
	}
	this->lambda = max(0.0,(E*sigma_sq_debt + E_debt*sigma_sq)/(sigma_sq*sigma_sq_debt*theta));
	this->w_assets = max(0.0,min(1.0-(0.2 + floor),E*sigma_sq_debt/(E*sigma_sq_debt + E_debt*sigma_sq)));
	this->folio_volume = max(0.0,min(this->getWealth(0.8)*mlimit, this->getWealth()*this->lambda));

		// pow(pow((1+lambda*(E*w_assets+(1-w_assets)*E_debt)-(theta*lambda*lambda*(sigma_sq*w_assets+(1-w_assets)*sigma_sq_debt))),(1.0-theta))/FFR,1.0/theta));
	}
double Node::assetSum(){
	double t=0;
	for (auto &asset: assets)
	{
		t += asset.second;
	}
	return t;
}



void Node::buyAssets(double amt, int mat){
	this->assets.push_back(make_pair(mat, amt));
	// cout<<"bought "<<assets[assets.size()-1].second<<" assets"<<endl;
}
// void Node::payIR(){
	
// }

/**
 * Default prediction (logit model trained on previous simulations)
 * Trained on: 300 simulations, 0 interest rate, balanced
 * LASSO regularization
 * Inverse of regularization strength = 0.0007
 * 
 * return -1.0 if believed to be called twice in one period
 */
double Node::defaultProb(){
	// coeff_deposits = 0.000982;
	// coeff_cash = -0.003140;
	// coeff_assets = 0.033246;
	// coeff_wealth = -0.128509;
	// coeff_leverage = 0.001108;
	// coeff_wealth_lag = -0.005691;
	// coeff_deposits_lag = 0.000325;
	// coeff_cash_lag = -0.001802;
	// coeff_assets_lag = 0.005801;
	
	double fitted = 0.000982 * this->deposits - 0.00314 * this->getCash() + 0.033246 * this->sumAssets()
			- 0.128509 * this->getWealth() + 0.001108 * this->getLeverage()
		    - 0.005691 * this->lag_wealth + 0.000325 * this->lag_deposits 
		    - 0.001802 * this->lag_cash * + 0.005801 * this->lag_sumAssets;

	double pred = exp(fitted) / (1 + exp(fitted)); // convert to probability

	// SAFER
	// check if lags and current differs in case it is called twice in one period
	// if (this->updateLags()){
	// 	cout<<"something is wrong with Node:defaultProb";
	// 	return -1.0;
	// }else{
	// 	return output;
	// }

	// LESS SAFE
	this->updateLags();

	cout << this->deposits << this->getCash() << this->sumAssets() << this->getWealth() << this->getLeverage()
		    << this->lag_wealth << this->lag_deposits << this->lag_cash << this->lag_sumAssets; // debug

	cout << this->nodeId << " Fitted=" << fitted << " DefaultPred=" << pred << endl; //debug

	return pred;
}

/**
 * update the lags if lags differ from current
 * I want to figure out a way to check whether the variables are updated 
 */
bool Node::updateLags(){
	cout << this->nodeId << " UpdateLags: " << "W=" << this->getWealth() << " D=" << this->deposits 
		<< " C=" << this->getCash() << " A=" << this->sumAssets() << endl; // debug

	// check if they are updated
	if (this->lag_wealth == this->getWealth() && this->lag_deposits == this->deposits
		&& this->lag_cash == this->getCash() && this->lag_sumAssets == this->sumAssets()){
		cout << this->nodeId << " NOT updated" << endl;
		return false; 
	}else{
		this->lag_wealth = this->getWealth();
		this->lag_deposits = this->deposits;
		this->lag_cash = this->getCash();
		this->lag_sumAssets = this->sumAssets();
		cout << this->nodeId << " YES updated" << endl;
		return true;
	}
}
