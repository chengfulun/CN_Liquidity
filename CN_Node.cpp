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


double Node::assetDemand(double price, double fundamental, double curr_price){
	double target = this->getWealth(curr_price) * (asset_aggression * )
	return this->asset_aggression * (price - fundamental) * this->getWealth(curr_price)
	- 
}

double Node::collateralValue(double CR, double price, bool useCredit){
	if (useCredit){
		double max_extrapolate = this->->getWealth(price) / CR;
		double investment_level = this->credit_target + this->asset_target;
		double max_payment = this->maxCredit(price) + this->getScrip();
		double result = 0;
		// change assetReturn to total investment return
		if (max_extrapolate < investment_level and investment_level < max_payment){
			result = this->assetReturn * (max_extrapolate - investment_level) / max_extrapolate;
		}
		if (max_extrapolate > investment_level and investment_level > max_payment){
			result = this->assetReturn * (max_extrapolate - investment_level) / max_extrapolate;
		}
		return result;
	}

	return CR * this->wealth_return;	
}

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

double Node::maxCredit(double haircut){
	double wealth = this->getWealth(haircut);
	double result = 0.0;
	double CR;
	vector <pair<double,double>> edges;
	for (auto &edge : atomicEdge_in){
		int tempix = edge.second->singleCreditIndex;
		double capacity = edge.second->capacity;
		CR = edge.second->originEdge->singleCreditEdges[tempix]->collateralRate;
		if (not edge.second->isDebt and wealth > 0){
			edges.push_back (std::make_pair(CR,capacity));
			// result += capacity;
			// wealth -= capacity * CR;
		}
	}
	sort(edges.begin(),edges.end());
	for (auto &e : edges){
		result += e.second;
		wealth -= e.second * e.first;
		if (wealth < 0){
			result -= wealth / e.first;
			break;
		}
	}
	return result;
}

double Node::targetCredit(){
	double credit_needed = this->credit_target + this->asset_target 
	- this->;
}


void Node::makeLeveraged(bool status, double haircut){
	if (status){
		this->leveraged = true;
		double ratio = this->getWealth(haircut) / (this->getDebt() + deposits);
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

double Node::expected_credit_revenue(){
	
}

double Node::getCollateral(double dRate){
	double temp = 0;
	double amt = 0;
	for (auto &aOut : atomicEdge_out){
		if (aOut->isDebt){
			temp += aOut->CR() * aOut->capacity;
			amt += aOut->capacity;
		}
	}

	// for (auto &eIn : edge_in){
	// 	for (auto &ceIn : eIn.second->singleCreditEdges){
	// 		double t = (ceIn->credit_max - ceIn->credit_remain->capacity);
	// 		amt += t;
	// 		temp += t*ceIn->collateralRate;
			// cout<<"debt: "<<t<<"at crate "<<ceIn->collateralRate<<endl;
	// 	}
	// }
	temp += deposits * dRate;
	// cout<<this->getNodeId()<<" total debt: "<<amt<<endl;
	// cout<<this->getNodeId()<<" total required collateral: "<<temp<<endl;
	return temp;
}

double Node::getScrip(){
	if(defaulted){
		return 0.0;
	}
	double temp = 0.0;
	for (auto &dIn : atomicEdge_in){
		if (dIn.second->isDebt){
			temp += dIn.second->capacity;
			// * (1.0 + dIn.second->interest_rate);
		}
	}
	return temp;
}

double Node::getCash(){
	if(defaulted){
		return 0.0;
	}
	double temp = 0.0;
	for (auto &dIn : atomicEdge_in){
		if (dIn.second->isDebt and dIn.second->nodeTo->isMarket){
			temp += dIn.second->capacity;
		}
	}
	// cout<<nodeNum<<endl;
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


double Node::getDemand(double price, double assetR){
	
}

double Node::getWealth(double haircut){
	if(defaulted){
		return 0.0;
	}
	// if(v>0){
	// // cout<<"node "<<nodeId<<" IOUs "<<this->getScrip()<< " debt "<<this->getDebt()<<" assets "<<this->sumAssets()<<" deposits "<<this->deposits<<" credit "<<getCredit()<<endl;		
	// }
	return this->getScrip() - this->getDebt() + this->sumAssets() * haircut - this->deposits;
	// double temp = 0;
	// cout<<nodeNum<<endl;
}

// double Node::getMaxPay(){
// 	double temp = 0.0;
// 	vector<pair<double,double>> credits;
// 	for (auto &eIn : atomicEdge_in){
// 		double cap = eIn.second->capacity;
// 		double cr = eIn.second->collateralRate
// 		if (eIn.second->isDebt){
// 			temp += ;
// 		}
// 		else{
// 			credits.push_back(make_pair())
// 		}
// 	}
// }



void Node::updateAssetReturn(){
	double result = 0.0;
	double sum_inv = 0.0;
	for (int i=0;i<asset_returns.size();++i){
		result+= asset_returns[i];
		sum_inv += asset_investments[i];
	}
	this->asset_return = result / sum_inv;
}

void Node::updateCreditReturn(){
	double result = 0.0;
	double sum_inv = 0.0;
	for (int i=0;i<credit_returns.size();++i){
		result+= credit_returns[i];
		sum_inv += credit_investments[i];
	}
	this->credit_return = result / sum_inv;
}

void Node::updateWealthReturn(){
	double result = 0.0;
	double sum_w = 0.0;
	for (int i=0;i<wealth_history.size();++i){
		result+= credit_returns[i] + asset_returns[i] - credit_investments[i] - credit_returns[i];
		sum_w += wealth_history[i];
	}
	this->wealth_return = result / sum_w;	
}

void Node::updateInvestmentReturn(){
	double result = 0.0;
	double sum_inv = 0.0;
	for (int i=0;i<asset_returns.size();++i){
		result += asset_returns[i] + credit_returns[i];
		sum_inv += credit_investments[i] + asset_investments[i];
	}
	this->investment_return = result / sum_inv;
}

void Node::updateInterestPaid(){
	double result = 0.0;
	double sum_inv = 0.0;
	for (int i=0;i<asset_leverage_interest.size();++i){
		result += asset_leverage_interest[i];
		sum_inv += asset_investments[i];
	}
	this->interest_paid = result / sum_inv;	
}

void Node::updateAllReturns(){
	this->updateAssetReturn();
	this->updateCreditReturn();
	this->updateWealthReturn();
	this->updateInvestmentReturn();
	this->updateInterestPaid();
}

void Node::updateAllHist(double price, double assetR, double creditR,double creditI,double assetI, double IR, int length){
	if (this->wealth_history.size() < length){
		wealth_history.erase(wealth_history.begin());
		asset_returns.erase(asset_returns.begin());
		credit_returns.erase(credit_returns.begin());
		asset_investments.erase(asset_investments.begin());
		credit_investments.erase(credit_investments.begin());
		asset_leverage_interest.erase(asset_leverage_interest.begin());
	}
	wealth_history.push_back(this->getWealth(price));
	asset_returns.push_back(assetR);
	credit_returns.push_back(creditR);
	asset_investments.push_back(assetI);
	credit_investments.push_back(creditI);
	asset_leverage_interest.push_back(IR);
}

double Node::getScrip(){
	if(defaulted){
		return 0.0;
	}
	double temp = 0.0;
	for (auto &dIn : atomicEdge_in){
		if (dIn.second->isDebt){
			temp += dIn.second->capacity;
			// * (1.0 + dIn.second->interest_rate);
		}
	}
	return temp;
}

double Node::creditReturn(double DR){
	if(defaulted){
		return 0.0;
	}
	double sumReturns = 0.0;
	for (auto &dIn : atomicEdge_in){
		if (dIn.second->isDebt){
			double valRate = dIn.second->interest_rate * (1.0 - DR)
			- (1.0 - dIn.second->CR()) * DR;
			sumReturns += dIn.second->capacity * valRate;
		}
	}
	return sumReturns;
}

double Node::credit2Return(double DR){
	if(defaulted){
		return 0.0;
	}
	double sumReturns = 0.0;
	for (auto &dIn : atomicEdge_in){
		if (dIn.second->isDebt){
			double valRate = dIn.second->interest_rate * (1.0 - DR)
			- (1.0 - dIn.second->CR()) * DR;
			sumReturns += dIn.second->capacity * valRate * valRate;
		}
	}
	return sumReturns;
}


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
	if (defaulted){
		return 0.0;
	}
	double temp = 0;
	for (auto &dOut : atomicEdge_out){
		if (dOut.second->isDebt && not dOut.second->nodeFrom->defaulted){
			temp += dOut.second->capacity;
			// * (1.0 + dOut.second->interest_rate);
		}
	}
	return temp;
}

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
	this->folio_volume = max(0.0,min(this->getWealth()*mlimit, this->getWealth()*this->lambda));
	
	this->asset_target = folio_volume * w_assets;
	this->credit_target = folio_volume * (1.0 - w_assets);
		// pow(pow((1+lambda*(E*w_assets+(1-w_assets)*E_debt)-(theta*lambda*lambda*(sigma_sq*w_assets+(1-w_assets)*sigma_sq_debt))),(1.0-theta))/FFR,1.0/theta));
	}

void Node::postCashUpdate(){
	double max_payment = this->maxCredit(1.0) + this->getScrip() - this->mReserve;
	this->folio_volume = min(this->folio_volume,max_payment);
	this->lambda = this->folio_volume / this->getWealth();
	this->asset_target = folio_volume * w_assets;
	this->credit_target = folio_volume * (1.0 - w_assets);
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