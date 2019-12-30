// CN_CreditNet.C 
#define NON_DEBUG_MODE

#include "CN_Constants.h"
#include "CN_CreditNet.h"
#include <string>
#include <queue>
#include <vector>
#include <iostream>
#include <math.h>
#include <cmath>

using namespace std;
// using namespace dlib;

extern CredNetConstants credNetConstants;

double phi(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
    return 0.5*(1.0 + sign*y);
}

double uniform_cdf(double x, double T){
	if (x <= -T){
		return 0.0;
	}
	if (x >= T){
		return 1.0;
	}
	return (x + T) / (2.0 * T);
}

CreditNet::CreditNet(double dVol, int finNumT, double precision, 
	int marketId, double initR,double initVol, 
	double deposit_rate, double deposit, double haircut, double mReserve, double mLimit,
	int valueBins, double EAR, double asset_vol, int defaulted_periods, double explore_boost, 
	bool basel, double max_lev, bool verb) : Graph(finNumT, marketId){
	this->precision = precision;
	this->transactionCounter = 0;
	this->marketId = marketId;
	this->initVol = initVol;
	this->initR = initR;
	this->deposit_rate = 0.0;
	this->haircut = haircut;
	// this->asset_maturity = maturity;
	this->valueBins = valueBins;
	this->verb = verb;
	this->mReserve = mReserve;
	this->expected_asset_return = EAR;
	this->expected_deposits = deposit;
	this->defaulted_periods = defaulted_periods;
	this->explore_boost = explore_boost;
	this->mLimit = mLimit;
	this->asset_volatility = asset_vol;

	this->log_mu = log(EAR*EAR/sqrt(EAR*EAR + asset_vol));
	// parameter sigma
	this->log_vol = sqrt(log(1+asset_vol/(EAR*EAR))); 	
	this->debt_coc_total = 0.0;

	this->dVol = dVol;
	this->dRevert = deposit_rate;

	this->basel = basel;
	this->max_lev = max_lev;
	// returns.assign(nodeNum, initR);
	// this->wealths.assign(nodeNum, 0.0);
	// credits.assign(nodeNum, 0.0);
	// credits_last.assign(nodeNum, 0.0);
	rDefaults.assign(nodeNum, 0.0);
	cDefaults.assign(nodeNum, 0.0);
	dDefaults.assign(nodeNum, 0.0);
	aDefaults.assign(nodeNum, 0.0);
	fDefaults.assign(nodeNum, 0.0);
	liquid_loss.assign(nodeNum, 0.0);
	default_wealth.assign(nodeNum, 0.0);
	deposit_losses.assign(nodeNum, 0.0);
	creditor_losses.assign(nodeNum, 0.0);
	credit_utils.assign(nodeNum,0.0);
	debt_coc.assign(nodeNum,0.0);
	asset_coc.assign(nodeNum,0.0);
	asset_amt.assign(nodeNum,0.0);
	debt_amt.assign(nodeNum,0.0);
	default_rates.assign(nodeNum, 0.0);
	credit_requests.assign(nodeNum,0.0);
	credit_premiums.assign(nodeNum,0.0);
	credit_utils.assign(nodeNum,0.0);
	debt_vol.assign(nodeNum,0.0);
	cReturns_real.assign(nodeNum,0.0);
	debt2_amt.assign(nodeNum,0.0);
	ntrials = 0;

	vector<int> initThis;
	for (int j=0;j<nodeNum - 1;j++){
		initThis.push_back(-1);
	}

	for (int i=0;i<nodeNum - 1;i++){		
		active_set_statuses.push_back(initThis);
	}
	// credits_sq.assign(nodeNum, 0.0);

	cReturns.assign(nodeNum, 0.0);
	c2Returns.assign(nodeNum,0.0);
	// cReturns_sq.assign(nodeNum, 0.0);
	// credits_cross.assign(nodeNum, 0.0);
	// volatilities.assign(nodeNum, initVol);
	for (int i=0;i<nodeNum;i++){
		this->nodes[i]->defaulted = false;
	}
}

CreditNet::CreditNet() : Graph(){}


vector<int> rix(int l){
	vector<int> v(l) ; // vector with 100 ints.
	iota (std::begin(v), std::end(v), 0); // Fill with 0, 1, ..., 99.
	random_shuffle ( v.begin(), v.end() );  // in place no extra array	
	return v;
}


void CreditNet::updateBasel(){
	double DR_all = total_dr / total_wealth;
	if (total_wealth == 0){
		DR_all = 0.0;
	}
	double sigmoid = 1 / (1+ exp(9.21 -614.51 * DR_all));
	this->mLimit = sigmoid + (1-sigmoid) * max_lev;
}


void CreditNet::record_values(){
	this->all_values.clear();
	for (auto &e : this->atomicEdges){
		if(not e.second->nodeTo->isMarket and not e.second->nodeTo->defaulted and not e.second->nodeFrom->defaulted){
			this->all_values.push_back(e.second->getSendValue(e.second->interest_rate));
			this->all_values.push_back(e.second->getReceiveValue(e.second->interest_rate));

		}
	}
}

void CreditNet::print_values(){
	for (auto &p : this->all_values){
	}
}

void CreditNet::init_lambdas(){
	for (int k = 0; k< nodeNum-1; k++){
		this->nodes[k]->getLambda(expected_asset_return - 1,asset_volatility * asset_volatility, initR,initVol,mLimit,1.0);
	}
}

void CreditNet::updateLambdas(){
	double AR;
	double CR;
	double CV;

	double all_AR = this->expected_asset_return -1 + asset_coc_all;
	double AV = this->asset_volatility;
	double all_CR = max(this->initR, creturn_all + debt_coc_all);
	double all_CV = max(this->initVol,creturn2_all);

	for (int k = 0; k< nodeNum-1; k++){
		AR = this->expected_asset_return -1 + asset_coc[k] / asset_amt[k];
		if (asset_amt[k] == 0){
			AR = expected_asset_return -1;
		}
		AR = max(AR,0.0);
		CR = cReturns[k]/debt_vol[k] + debt_coc_all;
			// debt_coc[k]) / debt_vol[k];
		if(debt_amt[k] == 0){
			CR = 0.0;
		}
		CR = max(this->initR,CR);

		CV = c2Returns[k] / debt2_amt[k];
		if (debt_amt[k] == 0){
			CV = initVol;
		}
		CV = max(this->initVol,CV);
		// multiplier limit is the leverage limit as a multiplier on wealth


		if(this->nodes[k]->new_flag > 0){
			AR = expected_asset_return - 1;
			CR = max(initR,all_CR);
			CV = max(all_CV,initVol);
		}
		this->nodes[k]->getLambda(AR, AV*AV, CR, CV, mLimit, 1.0);
		if(verb and AR != this->expected_asset_return - 1){
			cout<<"asset R: "<<AR<<" asset V: "<<AV<<" credit R: "<<CR<<" credit V: "<<CV<<endl;
		}
		if(verb){
			cout<<"node "<< k << " theta "<<this->nodes[k]->theta<<" portfolio size "<<this->nodes[k]->folio_volume<<" lambda "<<this->nodes[k]->lambda<<" assets "<<this->nodes[k]->w_assets<<" wealth "<<this->nodes[k]->getWealth(haircut)<<endl;
			cout<<"want assets "<<this->nodes[k]->getWealth(haircut)*this->nodes[k]->lambda*this->nodes[k]->w_assets<<endl;
		}
	}
}




void CreditNet::update_default_rates(){
	// get overall dr for basel
	this->total_dr = 0.0;

	for (int j = 0; j<nodeNum - 1; j++){
		this->default_rates[j] = this->defaultRate(j);
		this->total_dr += default_rates[j] * this->nodes[j]->getWealth(haircut);
		total_wealth += this->nodes[j]->getWealth(haircut);
		if(fabs(this->defaultRate(j)) > 1){
			cout<<"WARNING HIGH DR"<<endl;
		}
		if(verb){
			cout<<"nodenum: "<<j<<" default rate: "<<this->defaultRate(j)<<endl;			
		}

	}

	// total_dr = total_dr / dr_weight;


	for (auto &e : this->atomicEdges){
		double DR = this->default_rates[e.second->nodeTo->nodeId];
		e.second->updateDR(DR);
	}
}

double CreditNet::defaultRate(int nodeIx){
	if(nodeIx == marketId or this->nodes[nodeIx]->getWealth(haircut) > this->nodes[nodeIx]->sumAssets() or this->nodes[nodeIx]->sumAssets()==0){
		return 0.0;
	}
	double assetZero = 1.0 - this->nodes[nodeIx]->getWealth(haircut) / this->nodes[nodeIx]->sumAssets();

	if(assetZero <= 0){
		cout<<"node: "<<nodeIx<<" asset0: "<<assetZero<<" assets: "<<this->nodes[nodeIx]->sumAssets()<<" wealth: "<<this->nodes[nodeIx]->getWealth(haircut)<<endl;		
	}

	double insolventProb = phi((log(assetZero) - log_mu) / (sqrt(2.0) * log_vol));
	return insolventProb;
}


void CreditNet::updateCollateralValues(){
	double AR = this->expected_asset_return - 1 + asset_coc_all;
	double CR = creturn_all + debt_coc_all;
	for (auto &edge : this->atomicEdges){
		int k = edge.second->nodeTo->nodeId;
		AR = this->expected_asset_return -1 + asset_coc[k] / asset_amt[k];
		if (asset_amt[k] == 0){
			AR = expected_asset_return -1;
		}
		AR = max(AR,0.0);
		CR = cReturns[k]/debt_vol[k] + debt_coc_all;
			// debt_coc[k]) / debt_vol[k];
		if(debt_vol[k] == 0){
			CR = 0.0;
		}
		CR = max(this->initR,CR);

		edge.second->updateCollateralValue(CR,AR);
	}
}


void CreditNet::updateLineRequests(){
	for (int i = 0; i < nodeNum-1; i++){
		if (not this->nodes[i]->defaulted){
			double mult = min(mLimit,this->nodes[i]->lambda);
			double target = mult * this->nodes[i]->getWealth(haircut) + this->nodes[i]->reserveT;
			double payable = this->nodes[i]->getScrip() + this->nodes[i]->sumAssets() * expected_asset_return;
			this->credit_requests[i] = max(0.0,target - payable);
			if(verb){
				cout<<"node "<<i<<" requests "<<this->credit_requests[i]<<endl;
			}
		}
	}
}

void CreditNet::updateCreditUtil(){
	for (int i = 0; i<nodeNum-1; i++){
		credit_utils[i] = 0.0;
		if (not this->nodes[i]->defaulted){
			double issued = 0.0;
			double used = 0.0;
			for (auto &eOut : this->nodes[i]->edge_out){
				if (not eOut.second->nodeTo->isMarket){
					for(auto &it: eOut.second->singleCreditEdges){
						issued += it->credit_max;
						for(auto &d: it->debt_current){
							used += d.second->capacity;					
						}
					}
				}
			}
			credit_utils[i] = used / issued;
			if (issued == 0){
				credit_utils[i] = 0.0;
			}
		}

	}
}

void CreditNet::updateCreditPremiums(double thresh, double uprate, double downrate){
	updateCreditUtil();
	for (int i=0;i<nodeNum-1; i++){
		// updateNodePremium(thresh, uprate,downrate,credit_utils[i]);

		if (credit_utils[i] > thresh){
			this->nodes[i]->credit_premium *= uprate;
		}
		if (credit_utils[i] < thresh - 0.2){
			this->nodes[i]->credit_premium *= downrate;
		}
	}
}

// which nodes are in active set? can call while recording debt utilization
void CreditNet::updateCreditSets(double thresh, int inactive_periods){
	// each credit edge has a utilization, which should be used to update the nxn matrix of active statuses
	// active statuses are a n-vector for each node
	

	// double thresh = 0.25;
	// int inactive_periods = 2;

	for (auto &e : this->atomicEdges){
		if(e.second->nodeTo->isMarket){
			continue;
		}
		double UR = e.second->originEdge->singleCreditEdges[e.second->singleCreditIndex]->getUtil();
		int creditor = e.second->nodeFrom->nodeId;
		int borrower = e.second->nodeTo->nodeId;
		if (e.second->isDebt){
			continue;
		}
		// need node active status, not edge status
		if (active_set_statuses[creditor][borrower] == inactive_periods and UR < thresh){
			active_set_statuses[creditor][borrower] == -1;			
		}
		// if (active_set_statuses[creditor][borrower] == inactive_periods and UR >= thresh){
		// 	active_set_statuses[creditor][borrower] == 0;			
		// }
		if (active_set_statuses[creditor][borrower] < inactive_periods 
			and active_set_statuses[creditor][borrower] >= 0 
			and UR < thresh){
			active_set_statuses[creditor][borrower] += 1;			
		}
		// if (active_set_statuses[creditor][borrower] < inactive_periods 
		// 	and active_set_statuses[creditor][borrower] >= 0 
		// 	and UR >= thresh){
		// 	active_set_statuses[creditor][borrower] = 0;
		// }
		// if (active_set_statuses[creditor][borrower] == -1
		// 	and UR >= thresh){
		// 	active_set_statuses[creditor][borrower] = 0;			
		// }
		if(UR >= thresh){
			active_set_statuses[creditor][borrower] = 0;			
		}
	}
}

// cycle through requested lines and update capacities
void CreditNet::updateCreditLimits(){
	for (int j = 0; j < nodeNum - 1; j++){
		if(not this->nodes[j]->defaulted){
			for(auto& it : this->nodes[j]->atomicEdge_out){
				int fromId = it.second->nodeFrom->nodeId;
				if(it.second->isDebt && it.second->capacity>0 && not this->nodes[fromId]->defaulted){
					unordered_map<int,vector<double>>::const_iterator got = this->nodes[j]->creditPayOut.find (fromId);
					int tempix = it.second->singleCreditIndex;
					double tempCR = it.second->originEdge->singleCreditEdges[tempix]->collateralRate;
					double amt = it.second->capacity;
					if ( got == this->nodes[j]->creditPayOut.end() ){
						vector<double> new_elem;
						new_elem.push_back(tempCR);
						new_elem.push_back(amt*(1.0+it.second->interest_rate));
						new_elem.push_back(amt);
						this->nodes[j]->creditPayOut[fromId] = new_elem;
						this->nodes[j]->debt_out += it.second->capacity;
						this->nodes[fromId]->debt_in += it.second->capacity;					
					}
					else{
						// cout<<"before "<<this->nodes[j]->creditPayOut[fromId].second<<endl;
						this->nodes[j]->creditPayOut[fromId][1] += amt*(1.0+it.second->interest_rate);
						this->nodes[j]->creditPayOut[fromId][2] += amt;
						this->nodes[j]->debt_out += it.second->capacity;
						this->nodes[fromId]->debt_in += it.second->capacity;
						// cout<<"after "<<this->nodes[j]->creditPayOut[fromId].second<<endl;
					}
					total_ir += amt*it.second->interest_rate;
					total_cr += amt*tempCR;
					total_debt_res += amt;
					ir_collateralized += tempCR * it.second->interest_rate * amt;
					ir_free += (1-tempCR) * it.second->interest_rate * amt;
					debt_collateralized += tempCR * amt;
					debt_free += (1-tempCR)*amt;

					if(verb){
						cout<<"cpayout "<<this->nodes[j]->creditPayOut[fromId][1]<<" from "<<j<<" to "<<fromId<<endl;
					}

					//clear all debt after recording
					it.second->route(it.second->capacity,it.second->interest_rate,this,this->transactionCounter++);
				}
			}
		}

		// do i need been_defaulted logic to keep track of ressurrections?
		// the below just keeps track of empirical returns for nodes the just defaulted

		// if( this->nodes[j]->defaulted && not this->nodes[j]->been_defaulted){
		// 	for(auto& it : this->nodes[j]->atomicEdge_out){
		// 		int fromId = it.second->nodeFrom->nodeId;

		// 		if(it.second->isDebt && it.second->capacity>0){
		// 			cReturns[j] -= it.second->capacity;
		// 				// cout<<"after "<<this->nodes[j]->creditPayOut[fromId].second<<endl;
		// 		}
		// 			// cout<<"cpayout "<<this->nodes[j]->creditPayOut[fromId].second<<" from "<<j<<" to "<<fromId<<endl;

		// 		it.second->route(it.second->capacity,it.second->interest_rate,this,this->transactionCounter++);
		// 	}
		// 	this->nodes[j]->been_defaulted = true;

		// }
	}
	

	
	// need sum of active requests for each node's active set
	// credit investment target
	// record historic utilization on each credit line



	for (int j = 0; j < nodeNum -1; j++){
		// cout<<"leveraged? "<<overleveraged<<endl;
		// if(this->nodes[j]->defaulted){
		// 	cout<<j<<" already defaulted"<<endl;
		// 	continue;
		// }
		if (j != marketId && not this->nodes[j]->defaulted ){
			// cout<<j <<" not defaulted "<<this->nodes[j]->defaulted<<endl;&&not this->nodes[j]->leveraged
			double active_request_sum = 0.0;			
			vector<double> active_requests(this->nodeNum-1, 0.0);
			vector<int> inactive_ix;

			for (int k = 0; k < nodeNum - 1; k++){
				if(k==j){
					continue;
				}
				double crq = this->credit_requests[k];
				if (active_set_statuses[j][k] >= 0){
					active_requests[k] += crq;
					active_request_sum += crq;
				}
				else{
					inactive_ix.push_back(k);
				}
				if(verb){
					cout<<"active status: "<<active_set_statuses[j][k]<<" from node "<<j<<" to node "<<k<<endl;					
				}
			}

			// double w = wealths[j];
			// int overleveraged = checkCollateral(j);

			// }

			// double ir = (credNetConstants.uniformIntDistribution(
			// 	credNetConstants.globalGenerator) % numIR);
			// cout<<"random: "<<ir<<endl;	
				// this->addUnitEdge(nodes.find(i)->second, nodes.find(j)->second, ir, rand()%2);
			
			double inactive_budget = max(0.0,this->nodes[j]->credit_target - active_request_sum);
			double active_budget = this->nodes[j]->credit_target - inactive_budget;

			if(verb){
				cout<<"node "<<j<<"'s' inactive budget "<<inactive_budget<<" active budget "<<active_budget<<endl;				
			}
			vector<int> rand_ix = rix(nodeNum);
			for (int i = 0; i< nodeNum; i++){
				int ii = rand_ix[i];
			// cout<<cShares[i]<<endl;
				if(ii!= j && ii!= marketId and j!= marketId and not this->nodes[ii]->defaulted){
					double cap;					
					if(active_set_statuses[j][ii] >=0){
						// cShares[i] * this->nodes[j]->folio_volume*(1.0-this->nodes[j]->w_assets);
						// if(cap<0){cout<<cap<< " is credit"<<endl;}
						if (active_request_sum == 0.0){
							cap = active_requests[ii];
						}
						else{
							cap = active_requests[ii] * active_budget / (active_request_sum);
						}
						this->modifyCredit(j, ii, max(0.0,cap));
						active_budget -= cap;
					}
					else{
						// note that entire budget is not always spent
						cap = this->explore_boost * max(0.0,this->credit_requests[ii]);
						if(inactive_budget - cap > 0.0){
							this->modifyCredit(j, ii, cap);
							inactive_budget -= cap;					
						}
					}
					total_credit += cap;
					total_cr_extended += this->nodes[j]->CR_base * cap;
					// this->nodes[i]->print();
					// credits_sq[i] += cap*cap;

					// cout<<credits[i]<<endl;
				}
			}

			double extra_per_inactive = max(0.0,inactive_budget / (inactive_ix.size()) - 0.00000001);
			double extra_per_active = max(0.0, active_budget / (nodeNum - 1.00000001 - inactive_ix.size()));

			for(int i = 0;i<nodeNum - 1;i++){
				double cap;
				if(active_set_statuses[j][i] < 0 and inactive_budget > 0 and i!=j and not this->nodes[i]->defaulted){
					this->addCredit(j,i,extra_per_inactive);
					// cout<<"added "<<extra_per_inactive<<endl;
				}
				if(active_set_statuses[j][i] >= 0 and active_budget > 0 and i!=j and not this->nodes[i]->defaulted){
					this->addCredit(j,i,extra_per_active);
				}
			}
			
			


			// for(auto& it : this->nodes[j]->atomicEdge_in){
			// 	if(it.second->isDebt){
			// 		int toId = it.second->nodeTo->nodeId;
			// 		unordered_map<int,double>::const_iterator got = this->nodes[j]->creditPayIn.find (toId);

			// 		if ( got == this->nodes[j]->creditPayIn.end() )
			// 			this->nodes[j]->creditPayIn[toId] = it.second->capacity*(1+it.second->interest_rate);
			// 		else
			// 			this->nodes[j]->creditPayIn[toId] += it.second->capacity*(1+it.second->interest_rate);
			// 	}
			// }

		}

	}	
}

void CreditNet::updateCreditTerms(){
	for (auto &e : this->atomicEdges){
		if(e.second->isDebt or e.second->nodeFrom->isMarket or e.second->nodeTo->isMarket
			or e.second->nodeTo->defaulted){
			continue;
		}
		double DR = default_rates[e.second->nodeTo->nodeId];

		if (DR > 0.70){
			e.second->originEdge->singleCreditEdges[e.second->singleCreditIndex]->credit_max -= e.second->capacity;
			e.second->capacity = 0;
			continue;
		}
		double CR = e.second->nodeFrom->CR_base;
		double ir_set = min(0.2,(DR + e.second->nodeFrom->credit_premium - 
			CR * DR) / ((1-DR)));
		if(ir_set>=10 and verb){
			cout<<"IR ERROR"<<endl;
			cout<<"dr "<<DR<<" prem "<<e.second->nodeFrom->credit_premium<<endl;
			cout<<"node: "<<e.second->nodeFrom->nodeId<<" wants IR to be: "<<ir_set<<" for node "<<e.second->nodeTo->nodeId<<endl;			
		}
		e.second->originEdge->singleCreditEdges[e.second->singleCreditIndex]->collateralRate = CR;
		e.second->originEdge->singleCreditEdges[e.second->singleCreditIndex]->credit_interest_rate = ir_set;
		e.second->interest_rate = ir_set;
		total_ir_extended += ir_set * e.second->capacity;
	}
}


int CreditNet::makeInvest(bool forced, bool verbose){

	int leverages = 0;
	vector<int> randix  = rix(	nodeNum );

	for(int k = 0;k<nodeNum-1; k++){
		if (this->nodes[k]->defaulted){
			this->nodes[k]->default_counter += 1;
			if (this->nodes[k]->default_counter == this->defaulted_periods){
				this->restore_defaulted(k);
				if(verb){
					cout<<"restored node "<<k<<endl;
					this->nodes[k]->printBalance();
				}
			}
		}
		if(this->nodes[k]->new_flag > 0){
			this->nodes[k]->new_flag -= 1;
		}
		this->nodes[k]->creditPayOut.clear();
	}
	this->setReserves();

	for(int k = 0; k < nodeNum-1;k++){

		this->activateReserves(k);
		this->deactivateReserves(k,1);		
	}

	for (int kk = 0; kk< nodeNum; kk++){
		int k = randix[kk];
		if (not this->nodes[k]->defaulted and k != marketId){
			// cout<<"wealth before: "<<this->nodes[k]->getWealth(haircut)<<endl;		
			if(this->nodes[k]->debt_in != 0 or this->nodes[k]->debt_out != 0){
				cout<<"debt record error"<<endl;
				cout<<"debt_in: "<<	this->nodes[k]->debt_in<<" debt_out: "<<this->nodes[k]->debt_out<<endl;
			}
			double violate_lev = checkCollateral(k);
			leverages += violate_lev;
			cDefaults[k] += violate_lev;
			// cout<<"leveraged? "<<violate_lev<<endl;
			// if(this->nodes[j]->defaulted){
			// 	cout<<j<<" already defaulted"<<endl;
			// 	continue;
			// }

	// forced defaults
			// if(forced && not this->nodes[k]->defaulted && this->fDefault > 0){
			// 	this->nodes[k]->makeDefault();
			// 	if(this->verb){
			// 		cout<<"made "<<k<<" default "<<endl;				
			// 	}
			// 	this->fDefault -= 1;
			// 	fDefaults[k] += 1;
			// }
				// cout<<j <<" not defaulted "<<this->nodes[j]->defaulted<<endl;
				// double w = wealths[j];
			
			// cout<<"test"<<endl;
			// cout<<"deactivated reserve amount "<<this->nodes[k]->mReserved<<" available cash "<<this->nodes[k]->getCash()<<endl;
			// if(overleveraged == 0){
			// this->activateReserves(k);
			// this->deactivateReserves(k,1);	

			double fail;
			double debt_cost;
			double asset_cost;

			std::tie(fail,asset_cost,debt_cost) = payAsset(k, this->marketId, max(0.0,this->nodes[k]->asset_target), "MAX_FLOW_SRC", transactionCounter++, "ASSET",deposit_rate,haircut,0.0);
			if(verb){
				cout<<"paid assets "<<fail<<endl;
			}
			// this->nodes[i]->assets.push_back(std::make_pair(maturity,iAmount));

			if(fail > 0){
				this->nodes[k]->buyAssets(fail);
				// cout<<"asset shortage "<<this->nodes[k]->asset_target - fail<<endl;
				// cout<<"got assets "<<fail<<endl;
				if (fail > this->nodes[k]->asset_target){
					// cout<<"bought too many assets ERROR"<<endl;
				}
			}
			if(fail == -1){
				cout<<"CPLEX failure"<<endl;
				// print();
			}
			if(fail == 0){
				this->nodes[k]->buyAssets(this->nodes[k]->asset_target);
				// cout<<"got assets successfully "<<this->nodes[k]->asset_target<<endl;				
			}
			//remember to clear
			debt_coc[k] += debt_cost;
			// debt_coc_total += debt_cost;

			if(asset_cost< 0 and verb){
				cout<<"cost recorded: "<<asset_cost<<endl;
			}
			asset_coc[k]+= asset_cost;
			// activateReserves(k);		
			//this->wealths[k] = this->nodes[k]->getWealth(haircut);
			if(k != marketId){
				// cout<<"wealth after: "<<this->nodes[k]->getWealth(haircut)<<endl;			
			}
			// deactivateReserves(k);					
		}



	}

	for(int dd = 0; dd < nodeNum - 1; dd++){
		this->nodes[dd]->debt_in = 0.0;
		this->nodes[dd]->debt_out = 0.0;
		if((not this->nodes[dd]->defaulted) and this->nodes[dd]->getWealth(haircut) <= 0.0){
			if(verb){
				cout<<dd<<"after asset negative wealth default "<<this->nodes[dd]->getWealth(haircut)<<endl;
				this->nodes[dd]->printBalance();
			}
			this->makeDefault(dd);
			// defaultC+=1;
			aDefaults[dd] +=1;
			double dlost = this->nodes[dd]->deposits - this->nodes[dd]->getCash() - (0.5*this->nodes[dd]->sumAssets());
			liquid_loss[dd] += (0.5*this->nodes[dd]->sumAssets());
			deposit_losses[dd] += max(0.0,dlost);
			creditor_losses[dd] += min(this->nodes[dd]->getDebt()+dlost,this->nodes[dd]->getDebt());
			// unwind(dd);


		}	
	}
	// cout<<"finished asset purchases"<<endl;
	update_default_rates();
	if(basel){
		updateBasel();
	}
	// cout<<"got dr"<<endl;
	updateReturns();
	// cout<<"got returns"<<endl;
	updateLambdas();
	// cout<<"got lambdas"<<endl;




	updateLineRequests();
	// cout<<"got lines"<<endl;
	updateCreditPremiums(0.95, 1.05, 0.9);
	// cout<<"got prem"<<endl;
	updateCreditSets(0.25, 2);
	// cout<<"got csets"<<endl;
	updateCreditLimits();
	// cout<<"got climits"<<endl;
	updateCreditTerms();
	// cout<<"got terms"<<endl;
	updateCollateralValues();

	record_values();
	credNetConstants.setValues(all_values,valueBins);

	// cout<<"cleared hist"<<endl;	
	// each credit edge has a utilization, which should be used to update the nxn matrix of active statuses
	// active statuses are a n-vector for each node
	
	ntrials += 1;

	total_defaults_sum += total_defaults;
	total_ir_sum += total_ir;
	total_cr_sum += total_cr;
	total_debt_res_sum += total_debt_res;
	total_wealth_sum += total_wealth;
	capital_loss_sum += capital_loss;
	total_dr_sum += total_dr;

	total_credit_sum += total_credit;
	total_cr_extended_sum += total_cr_extended;
	total_ir_extended_sum +=total_ir_extended;
	ir_collateralized_sum += ir_collateralized;
	ir_free_sum += ir_free;
	debt_collateralized_sum += debt_collateralized;
	debt_free_sum += debt_free;

	mLimit_sum += mLimit * total_wealth;
	return leverages;



	// for(int i=0;i<returns.size() ; i++){
	// 	cout<<returns[i]<<endl;
	// // }

	// credits_last=credits;

	// record unmodified credit, reset actual credit
	
}

void CreditNet::clear_history(){
	std::fill(debt_coc.begin(), debt_coc.end(), 0.0);
	std::fill(asset_coc.begin(), asset_coc.end(), 0.0);
	std::fill(asset_amt.begin(), asset_amt.end(), 0.0);
	std::fill(debt_amt.begin(), debt_amt.end(), 0.0);
	std::fill(cReturns.begin(), cReturns.end(), 0.0);
	std::fill(c2Returns.begin(), c2Returns.end(), 0.0);
	std::fill(cReturns_real.begin(), cReturns_real.end(), 0.0);
	std::fill(debt2_amt.begin(),debt2_amt.end(),0.0);
	std::fill(debt_vol.begin(),debt_vol.end(),0.0);	
	debt_coc_total = 0.0;
	total_ir = 0;
	total_cr = 0;
	total_dr = 0;
	total_debt_res = 0;
	total_wealth = 0;
	capital_loss = 0;
	total_defaults = 0;

	total_credit = 0;
	total_cr_extended = 0;
	total_ir_extended = 0;
	ir_collateralized = 0;
	ir_free = 0;
	debt_collateralized = 0;
	debt_free = 0;

	defaultable_periods = nodeNum-1;
}

void CreditNet::restore_defaulted(int k){
	this->nodes[k]->deposits = this->expected_deposits;
	this->nodes[k]->defaulted = false;
	this->nodes[k]->default_counter = 0;
	for (auto &eIn: this->nodes[k]->atomicEdge_in){
		if (1){
			eIn.second->capacity = 0;
		}
	}

	for (auto &eOut: this->nodes[k]->atomicEdge_out){
		if (eOut.second->isDebt){
			eOut.second->capacity = 0;
		}
		else if (not eOut.second->nodeTo->isMarket){
			eOut.second->capacity = 0;			
		}
	}

	grabCollateral(k,this->expected_deposits + this->init_wealth);
	this->nodes[k]->new_flag= 2;
	// cout<<"node restored"<<endl;
}

void CreditNet::resultsOut(){
	// cout<<"printing"<<endl;
	for (int i = 0; i< nodeNum - 1; i++){
		std::vector<double> adjrec(nodeNum - 1, 0.0);
		for (auto &eIn : this->nodes[i]->edge_in){
			if(not eIn.second->nodeFrom->defaulted){
				for (auto &cEdge : eIn.second->singleCreditEdges){
					for (auto &debtEdge : cEdge->debt_current){
						adjrec[eIn.second->nodeFrom->nodeId] += debtEdge.second->capacity;
					}
				}
			}
		}
		for (int j = 0; j < nodeNum - 1; j++){
			cout<< adjrec[j] << " ";
		}
		cout<<this->nodes[i]->theta<< " "<<rDefaults[i]
		<<" "<<aDefaults[i]<<" "<<dDefaults[i]<<" "<<cDefaults[i]<<" "<<this->nodes[i]->getWealth(haircut)
		<<" "<<this->nodes[i]->getDebt()<<" "<<
		this->nodes[i]->getCredit()<<" "<<
		this->nodes[i]->getScrip() - this->nodes[i]->getCash()<<" "<<
		this->nodes[i]->folio_volume*(1-this->nodes[i]->w_assets)<<" "<<this->nodes[i]->deposits<<" "<<
		this->nodes[i]->getCash()<< " "<<this->nodes[i]->sumAssets()<<endl;
	}
}



void CreditNet::resultsOut_1(){
	// cout<<"printing"<<endl;
	for (int i = 0; i< nodeNum - 1; i++){
		cout<<i<<" "<<ntrials<<" "<<this->nodes[i]->theta<< " "<<this->nodes[i]->CR_base<<" "<<rDefaults[i]
		<<" "<<aDefaults[i]<<" "<<dDefaults[i]<<" "<<cDefaults[i]<<" "<<this->nodes[i]->getWealth(haircut)
		<<" "<<this->nodes[i]->debt_out<<" "<<
		this->nodes[i]->getCredit()<<" "<<
		this->nodes[i]->getScrip() - this->nodes[i]->getCash()<<" "<<
		this->nodes[i]->credit_target<< " "<< this->nodes[i]->deposits<<" "<<this->nodes[i]->getCash()<< " "<<this->nodes[i]->sumAssets()<<
		" "<<this->nodes[i]->mReserved<<" "<<this->alpha<<" "<<liquid_loss[i]<<" "<<
		default_wealth[i]<< " "<<deposit_losses[i]<<" "<<
		this->default_rates[i]<<" "<<this->debt_coc_all<<" "<<this->asset_coc[i] / asset_amt[i]<<" "<<
		this->cReturns[i] / debt_vol[i]<<" "<<this->cReturns_real[i] / debt_vol[i]<<" "<<
		this->creditor_losses[i]<<this->expected_asset_return -1 + asset_coc[i] / asset_amt[i]<<" "<<
		cReturns[i]/debt_vol[i] + debt_coc_all<<" "<<
		c2Returns[i] / debt2_amt[i]<<" "<<this->default_rates[i]<<endl;
	}
}


void CreditNet::results_macro(){
	cout<<ntrials<<" "<<total_defaults<<" "<<
	defaultable_periods<<" "<<total_ir / total_debt_res<<" "<<
	total_debt_res<<" "<<total_cr / total_debt_res<<" "<<total_dr / total_wealth<<
	" "<<total_wealth<<" "<<capital_loss<<" "<<
	total_cr_extended / total_credit<<" "<<
	total_ir_extended / total_credit<<" "<<
	total_credit<<" "<<
	ir_collateralized / debt_collateralized<<" "<<
	ir_free / debt_free<<" "<<
	debt_collateralized<<" "<<debt_free<<" "<<mLimit
	<<endl;	// asset_coc_all<<" "<<debt_coc_all<<endl;
}

void CreditNet::results_macro_all(){
	cout<<-1<<" "<<total_defaults_sum<<" "<<
	defaultable_periods_sum<<" "<<total_ir_sum / total_debt_res_sum<<" "<<
	total_debt_res_sum<<" "<<total_cr_sum / total_debt_res_sum<<" "<<total_dr_sum / total_wealth_sum<<
	" "<<total_wealth_sum<<" "<<capital_loss_sum<<" "<<
	total_cr_extended_sum / total_credit_sum<<" "<<
	total_ir_extended_sum / total_credit_sum<<" "<<
	total_credit_sum<<" "<<
	ir_collateralized_sum / debt_collateralized_sum<<" "<<
	ir_free_sum / debt_free_sum<<" "<<
	debt_collateralized_sum<<" "<<debt_free_sum<<" "<<mLimit_sum / total_wealth_sum
	<<endl;
	// asset_coc_all<<" "<<debt_coc_all<<endl;
}

int CreditNet::shockPay(double alpha, bool crisis){
	this->alpha = alpha;
	// fill(cReturns.begin(), cReturns.end(), 0.0);
	// cReturns[fid] = 0.0;	
	// if(this->nodes[fid]->defaulted){
	// 	cout<<fid<<" already defaulted"<<endl;
	// 	return 0;
	// }
	// cout<<fid <<" not defaulted shock "<<this->nodes[fid]->defaulted<<endl;
	
	//deposit shock
	clear_history();

	for (int fid = 0; fid< nodeNum-1;fid++){
		if (this->nodes[fid]->defaulted){
			defaultable_periods -= 1;
		}
	}

	defaultable_periods_sum += defaultable_periods;

	vector<int>randix  = rix(	nodeNum );

	for (int fid1 = 0;fid1<nodeNum;fid1++){
		int fid = randix[fid1];
		if(fid!=marketId && not this->nodes[fid]->defaulted){
			double dShock = (( (credNetConstants.uniformDoubleDistribution(
								credNetConstants.globalGenerator)) - 0.5) * 2.0 ) * this->nodes[fid]->deposits;
			
			// double dNorm = (credNetConstants.normalDistribution(
			// 						credNetConstants.globalGenerator))*dVol + log_mu;
			

			// cout<<"dShock is "<<dShock<<endl;
			activateReserves(fid);
			if(not crisis || (crisis && this->cCount <= 0)){
				// double dShock = (credNetConstants.normalDistribution(
				// 								credNetConstants.globalGenerator)) * deposit_shock * this->nodes[fid]->deposits;
				// double dShock = ( (credNetConstants.uniformDoubleDistribution(
				// 					credNetConstants.globalGenerator)) - 0.5) * 2.0* deposit_shock * this->nodes[fid]->getWealth(1.0);
				
				this->nodes[fid]->deposits += dShock;
				// if(this->nodes[fid]->deposits<0){
				// 	dShock = dShock - this->nodes[fid]->deposits;
				// 	this->nodes[fid]->deposits = 0.0;
				// }
				if(verb){
					cout<<"node "<<fid<<" deposit shock "<<dShock<<" cash "<<this->nodes[fid]->getCash()<< " credit "<<this->nodes[fid]->getCredit()<<" assets "<< this->nodes[fid]->sumAssets()<<" leveraged "<<this->nodes[fid]->leveraged<<endl;
				}
				if(dShock >= 0){
					grabCollateral(fid,dShock);
				}
				else{
					double failedDeposit;
					double debt_cost;
					double asset_cost;
					if(verb){
						cout<<"pay deposit"<<endl;						
					}
					std::tie(failedDeposit,asset_cost,debt_cost) = payAsset(fid, this->marketId, -dShock, "MAX_FLOW_SRC", transactionCounter++, "DEBT",deposit_rate,haircut,0.0);					
					if (failedDeposit > 0.0 or failedDeposit == -1){
						if(failedDeposit == -1){
							failedDeposit = 0;
						}
						// cout<<"was defaulted "<< this->nodes[fid]->defaulted<<endl;
						double status = makeLiquidate(fid,-dShock - failedDeposit);
						liquid_loss[fid] += status;

						if(status >0){
							if(verb){
								cout<<fid << " couldn't pay " << -dShock - failedDeposit<< " is defaulted "<<this->nodes[fid]->defaulted<<endl;								
							}
							// this->nodes[fid]->makeDefault();
							// defaultC+=1;
							dDefaults[fid] +=1;
							double dlost = this->nodes[fid]->deposits - this->nodes[fid]->getCash() +(-dShock - failedDeposit) -status;
							default_wealth[fid] += max(0.0,this->nodes[fid]->getWealth(1.0));
							deposit_losses[fid] += max(0.0,dlost);
							creditor_losses[fid] += min(this->nodes[fid]->getDebt()+dlost,this->nodes[fid]->getDebt());
						}

						
					}
					
					// debt_coc_total += debt_cost;
					debt_coc[fid] += debt_cost;
				}				
			}
			// if (crisis && this->cCount>0){
			// 	if(verb){
			// 		cout<<"crisis at "<<fid<<endl;
			// 	}
			// 	if(verb){
			// 		cout<<"node "<<fid<<" deposit shock "<<this->nodes[fid]->deposits<<" cash "<<this->nodes[fid]->getCash()<< " credit "<<this->nodes[fid]->getCredit()<<" assets "<< this->nodes[fid]->sumAssets()<<endl;
			// 	}				
			// 	this->cCount += -1;
			// 	fDefaults[fid] +=1;										
			// 	double failedDeposit = payAsset(fid, this->marketId, this->nodes[fid]->deposits, "MAX_FLOW", transactionCounter++, "DEBT",deposit_rate, deposit_rate, haircut);
			// 	if (failedDeposit > 0){
			// 		// cout<<"was defaulted "<< this->nodes[fid]->defaulted<<endl;
			// 		double status = makeLiquidate(fid,this->nodes[fid]->deposits - failedDeposit);
			// 		liquid_loss[fid] += status;
			// 		if(status >0){
			// 			if(verb){
			// 				cout<<fid << " couldn't pay " << this->nodes[fid]->deposits - failedDeposit<< " is defaulted "<<this->nodes[fid]->defaulted<<endl;								
			// 			}
			// 			// this->nodes[fid]->makeDefault();
			// 			defaultC+=1;
			// 			dDefaults[fid] +=1;	
			// 			double dlost = this->nodes[fid]->deposits - this->nodes[fid]->getCash() +(this->nodes[fid]->deposits - failedDeposit) -status;
			// 			default_wealth[fid] += max(0.0,this->nodes[fid]->getWealth(1.0));
			// 			deposit_losses[fid] += max(0.0,dlost);
			// 			creditor_losses[fid] += min(this->nodes[fid]->getDebt()+dlost,this->nodes[fid]->getDebt());											

			// 		}
			
			// 	}
			// 	this->nodes[fid]->deposits = 0;				
			// }
			this->deactivateReserves(fid);
		}
	}



	//asset shock
	int defaultC = 0;
	for (int fid = 0;fid<nodeNum - 1;fid++){
		if(fid!=marketId && not this->nodes[fid]->defaulted){
			grabCollateral(fid, this->nodes[fid]->assets*(alpha));
			this->nodes[fid]->assets = 0;
		}
	}
			// auto i = begin(this->nodes[fid]->assets);
			// while( i != end(this->nodes[fid]->assets)){
			// 	int ix = i - this->nodes[fid]->assets.begin();
			// 	cout<<"asset ix "<<ix<<endl;
			// 	this->nodes[fid]->assets[ix].first -= 1;
			// 	cout<<"node "<< fid<<" maturity "<< this->nodes[fid]->assets[ix].first<<" amt "<<this->nodes[fid]->assets[ix].second<<endl;
			// 	if ( this->nodes[fid]->assets[ix].first == 0){
			// 		double rr = (credNetConstants.uniformIntDistribution(
			// 			credNetConstants.globalGenerator) % 7);
			// 		double alpha = (rr - 3)/(2/asset_volatility) + expected_asset_return;
			// 		cout<<"alpha is "<<alpha<<endl;
			// 		grabCollateral(this->nodes[fid]->nodeId, this->nodes[fid]->assets[ix].second*(1+alpha));
			// 		cout<<"asset paid off "<<this->nodes[fid]->assets[ix].second*(1+alpha)<<endl;
			// 		i = this->nodes[fid]->assets.erase(i);
			// 	}
			// 	else{
			// 		++i;
			// 	}
			// }


	// debt repayment

	randix = rix(nodeNum);

	for (int f = 0;f<nodeNum ;f++){
		int ff = randix[f];
		if( ff != marketId && not this->nodes[ff]->defaulted){
			this->activateReserves(ff);
			for (auto& cPay : this->nodes[ff]->creditPayOut){
				if(verb){
					cout<<"pay credit"<<endl;
					cout<<"c to pay "<<cPay.second[1]<<" from "<<ff<<" to "<<cPay.first<<endl;
					cout<<"cash: "<<this->nodes[ff]->getCash()<<endl;					
					this->nodes[ff]->printBalance();
				}
				double current_crate = cPay.second[0];
				double paymentFail;
				double debt_cost;
				double asset_cost;
				std::tie(paymentFail,asset_cost,debt_cost) = payAsset(ff, cPay.first, cPay.second[1], "MAX_FLOW_SRC", transactionCounter++, "DEBT",deposit_rate,haircut,current_crate);					

				if(paymentFail > 0){

					double left = cPay.second[1] - paymentFail;
					double status = payCash(ff, left);
					if(status > 0){
						capital_loss += cPay.second[1];
						if(verb){
							cout<<"node "<<ff<<" couldn't pay credit "<<cPay.second[1] - paymentFail<<" defaulted"<<endl;
							this->nodes[ff]->printBalance();								
						}						
						makeDefault(ff);
						default_wealth[ff] += max(0.0,this->nodes[ff]->getWealth(1.0));
						// 
						// double status = makeLiquidate(ff,cPay.second.second - paymentFail);
						// if(status >0){
						// 	// this->nodes[fid]->makeDefault();


						// 	defaultC+=1;
						rDefaults[ff] +=1;						

						// 	// grabCollateral(cPay.first,status);
						cReturns_real[ff] += paymentFail;											
						// }
						// else{
						// 	grabCollateral(cPay.first,cPay.second.second - paymentFail);
						// }						
						this->nodes[cPay.first]->debt_in -= cPay.second[2];
						double status = payCash(cPay.first,paymentFail);
					}
					else{
						grabCollateral(cPay.first,left);
						cReturns_real[ff] += cPay.second[1];
						this->nodes[ff]->debt_out -= cPay.second[2];
						// -= cPay.second.second;
						this->nodes[cPay.first]->debt_in -= cPay.second[2];
						}

				}
				else{
					// grabCollateral(cPay.first,cPay.second[1]);
					cReturns_real[ff] += cPay.second[1];	
					this->nodes[ff]->debt_out -= cPay.second[2];
					// -= cPay.second.second;
					this->nodes[cPay.first]->debt_in -= cPay.second[2];
				}
				cPay.second.clear();

				// cReturns_sq[fid] += cPay.second.second*cPay.second.second;
				// credits_cross[fid] += cPay.second.second * this->cap;
				if(not this->nodes[ff]->defaulted){

				}

				// -= cPay.second.second;					cPay.second.second = 0;					
				// cout<<"cpaid total "<<cReturns[ff]<<endl;
				// debt_coc_total += debt_cost;
				debt_coc[ff] += debt_cost;				
			}
			// if(this->nodes[ff]->defaulted){
			// 	this->nodes[ff]->been_defaulted += 1;
			// 	credits_last[ff]+= credits[ff];
			// 	cReturns[ff] += credits[ff] - ;
			// 	// returns[ff] = (cReturns[ff] / credits_last[ff]) - 1.0;
			// }
			this->deactivateReserves(ff);
		}		
	}



	// handle defaults
	for(int dd = 0; dd < nodeNum - 1; dd++){
		this->nodes[dd]->debt_in = 0.0;
		this->nodes[dd]->debt_out = 0.0;		
		if((not this->nodes[dd]->defaulted) and this->nodes[dd]->getWealth(haircut) <= 0.0){
			if(verb){
				cout<<dd<<" negative wealth default "<<this->nodes[dd]->getWealth(haircut)<<endl;
				this->nodes[dd]->printBalance();										
			}
			this->makeDefault(dd);
			// defaultC+=1;
			aDefaults[dd] +=1;
			double dlost = this->nodes[dd]->deposits - this->nodes[dd]->getCash() - (0.5*this->nodes[dd]->sumAssets());
			liquid_loss[dd] += (0.5*this->nodes[dd]->sumAssets());
			deposit_losses[dd] += max(0.0,dlost);
			creditor_losses[dd] += min(this->nodes[dd]->getDebt()+dlost,this->nodes[dd]->getDebt());
			// unwind(dd);


		}	
	}


	// modify investment targets to account for mReserve
	this->postCashUpdateLambda();
	setReserves();

	if(verb){
		cout<<"loop done"<<endl;
	}

	return defaultC;		
}

double CreditNet::payCash(int fid, double amt){
	double leftover = amt;
	for (auto&it : this->nodes[fid]->atomicEdge_in){
		if(it.second->nodeTo->isMarket and it.second->isDebt){
			double routable = max(0.0,min(it.second->capacity,leftover));	
			it.second->route(routable,0.0,this,this->transactionCounter++);
			leftover -= routable;
		}
		if (leftover<=0){
			it.second->originEdge->singleCreditEdges[it.second->singleCreditIndex]->credit_remain->route(-leftover,0,this,this->transactionCounter++);
			return 0;
		}
	}
	return leftover;
}


void CreditNet::makeDefault(int fid){
	if (verb){
		cout<<"defaulted node: "<<fid<<" epoch "<<ntrials<<endl;	
		cout<<"assets: "<<this->nodes[fid]->sumAssets()<<" shock "<<alpha<<endl;
		cout<<"cash: "<<this->nodes[fid]->getCash()<<endl;
		cout<<"debt: "<<this->nodes[fid]->getDebt()<<endl;
		cout<<"IOUs: "<<this->nodes[fid]->getScrip() - this->nodes[fid]->getCash()<<endl;
		cout<<"deposits: "<<this->nodes[fid]->deposits<<endl;
		cout<<"unpaid debt out: "<<this->nodes[fid]->debt_out<<endl;
		cout<<"unpaid debt in: "<<this->nodes[fid]->debt_in<<endl;		
	}
	total_defaults += 1;
	// this->nodes[fid]->print();
	this->nodes[fid]->defaulted = true;
	this->nodes[fid]->reserveT = 0.0;
	this->nodes[fid]->debt_reserve = 0.0;
	this->nodes[fid]->mReserved = 0.0;
	for (auto &eIn : this->nodes[fid]->atomicEdge_in){
		if(eIn.second->isDebt){
			// cout<<"pay defaulted node"<<endl;
			double owed = eIn.second->capacity;
			if(owed == 0){
				continue;
			}
			double cr_owed = eIn.second->CR();
			double paymentFail;
			double debt_cost;
			double asset_cost;

			std::tie(paymentFail,asset_cost,debt_cost) = payAsset(eIn.second->nodeFrom->nodeId, this->marketId, owed, "MAX_FLOW_SRC", transactionCounter++, "DEBT",deposit_rate,haircut,cr_owed);			
			if(paymentFail > 0){
				double left = payCash(eIn.second->nodeFrom->nodeId,owed-paymentFail);
				if (left > 0.5 * eIn.second->nodeFrom->assets){
					makeDefault(eIn.second->nodeFrom->nodeId);
					if(verb){
						// cout<<"node "<<fid<<" couldn't pay credit "<<owed - paymentFail<<" defaulted"<<endl;	
					}
					// defaultC+=1;
					rDefaults[fid] +=1;				
				}
				else{
					eIn.second->nodeFrom->assets -= 2 * (owed-paymentFail);
				}
				
				

					// grabCollateral(cPay.first,status);
				
				// else{
				// 	grabCollateral(cPay.first,cPay.second.second - paymentFail);
				// }
			}
			debt_coc[fid] += debt_cost;
		}
		eIn.second->capacity = 0.0;

	}
	this->nodes[fid]->debt_in = 0.0;
	this->nodes[fid]->debt_out = 0.0;

	for(auto &eOut : this->nodes[fid]->atomicEdge_out){
		if(eOut.second->isDebt){
			double owe = eOut.second->capacity * eOut.second->CR();
			capital_loss += eOut.second->capacity * (1-eOut.second->CR());
			grabCollateral(eOut.second->nodeFrom->nodeId,owe);
			this->cReturns_real[eOut.second->nodeFrom->nodeId] += owe;
		}
		if(not eOut.second->nodeTo->isMarket){
			eOut.second->capacity = 0.0;			
		}
	}
}

void CreditNet::setReserves(){
for (int k = 0; k< nodeNum - 1; k++){
	if(this->nodes[k]->defaulted){
		continue;
	}
	this->nodes[k]->reserveT = min(this->nodes[k]->deposits * mReserve,this->nodes[k]->getCash());
}
}


void CreditNet::postCashUpdateLambda(){
	for (int k = 0; k< nodeNum - 1; k++){
		if(this->nodes[k]->defaulted){
			continue;
		}
		this->nodes[k]->postCashUpdate(haircut);
	}
}

void CreditNet::updateReturns(){
	if(ntrials <2){
		for(int i = 0; i<nodeNum -1;i++){
			if(not this->nodes[i]->defaulted){
				cReturns[i] = this->initR;
				c2Returns[i] = this->initVol;
				debt_amt[i] = 1.0;
				debt2_amt[i] = 1.0;
				asset_amt[i] = 1.0;
				asset_coc[i] = 0.0;
				debt_coc[i] = 0.0;
				debt_vol[i] = 0.0;
			}
		}
		asset_coc_all = 0.0;
		debt_coc_all = 0.0;
		creturn_all = initR;
		creturn2_all = initVol;
	}
	else{
		double total_asset_coc = 0.0;
		double total_debt_coc = 0.0;
		double total_debt_value = 0.0;
		// amount of debt, so we can assign cash reserved for debt
		double total_debt_vol = 0.0;
		double total_debt2 = 0.0;

		// total cash assigned to funding debt
		total_debt = 0.0;
		// total_debt2 = 0.0;
		total_assets = 0.0;
		asset_coc_all = 0.0;
		debt_coc_all = 0.0;
		creturn_all = 0.0;
		creturn2_all = 0.0;

		for (int i = 0; i< nodeNum -1; i++){
			if(this->nodes[i]->defaulted){
				asset_amt[i] = 0;
				debt_amt[i] = 0;
				cReturns[i] = 0;
				c2Returns[i] = 0;
				debt_vol[i] = 0;
				continue;
			}
			asset_amt[i] = this->nodes[i]->sumAssets();
			asset_coc_all += asset_coc[i];
			total_assets += asset_amt[i];
			
			debt_amt[i] = this->nodes[i]->getScrip() - this->nodes[i]->getCash();
			debt2_amt[i] = this->nodes[i]->scrip2Sum();
			total_debt2 += debt2_amt[i];

			debt_coc_all += debt_coc[i];
			total_debt += debt_amt[i];
			debt_vol[i] = max(this->nodes[i]->debt_reserve,this->nodes[i]->getDebtCaps());
			total_debt_vol += debt_vol[i];

			cReturns[i] = this->nodes[i]->creditReturn(this->default_rates[i]);
			c2Returns[i] = this->nodes[i]->credit2Return(this->default_rates[i],cReturns[i] / debt_amt[i]);

			// double creturn = this->nodes[i]->creditReturn(this->default_rates[i]);
			// double c2return = this->nodes[i]->credit2Return(this->default_rates[i]);
			creturn_all += cReturns[i];
			creturn2_all += c2Returns[i];
			// cReturns[i] = creturn / debt_amt[i];
			// c2Returns[i] = c2return / debt_amt[i];


			// double ir_rate = 0;
			// double 

			// double delta;
			// double delta2;
			// // if(this->nodes[i]->defaulted){
			// // 	returns[i] = 0.0;
			// // 	volatilities[i] = 0.0;
			// // }
			// // else{
			// if(credits_last[i]<0.00001){
			// 	delta = 0.0;
			// }
			// // else if(aDefaults[i] == 1){
			// // 	delta = -returns[i];
			// // }
			// else{
			// 	delta = (cReturns[i] - credits_last[i])/(credits_last[i] + 0.00000000000001) - returns[i];

			// }
			// // if(cReturns[i] > 0){
			// 	// cout<<"credit returns "<<cReturns[i]<<" credit expenditure "<<credits_last[i]<<endl;
			// // }
			// if(cReturns[i]>credits_last[i]*2.0 && verb){
			// 	cout<<"credits violated "<<i<<endl;
			// }
			// // cout<<"leveraged "<<this->nodes[i]->leveraged<<" defaulted "<<this->nodes[i]->defaulted<<endl;
			// returns[i] = returns[i] + delta/(double) ntrials;
			// // returns[i] = cReturns[i]/credits[i];
			// if(verb){
			// 	cout<<"node "<<i<<" returns delta "<<delta<<" returns "<<cReturns[i]<<" credits "<<credits_last[i]<<" defaulted "<<this->nodes[i]->defaulted<<endl;
			// }
			// // if(this->nodes[i]->defaulted){
			// // 	delta2 = -
			// // }
			// delta2 = (cReturns[i] - credits_last[i])/(credits_last[i] + 0.00000000000001) - returns[i];
			// // cout<<"vols delta "<<delta2<<endl;

			// volatilities[i] = (volatilities[i]*(double) (ntrials-2) + (delta*delta2))/(double) (ntrials - 1);

			// }
		}
		asset_coc_all = asset_coc_all / total_assets;
		if(total_assets == 0){
			asset_coc_all = 0.0;
		}

		debt_coc_all = debt_coc_total / total_debt_vol;
		// total_debt_vol;
		if(total_debt_vol == 0){
			debt_coc_all = 0.0;
		}

		creturn_all = creturn_all / total_debt_vol;
		// total_debt_vol;
		if(total_debt_vol == 0){
			creturn_all = 0.0;
		}
		creturn2_all = creturn2_all / total_debt2;
		// total_debt_vol;
		if(total_debt2 == 0){
			creturn2_all = 0.0;
		}
	}
	// ntrials += 1;
}

double CreditNet::makeLiquidate(int fid, double amt){
	if(verb){
		cout<<"liquidating "<<amt<<" for node "<<fid<<endl;
	}
	double liquidate = amt * 2.0;
	if (liquidate > this->nodes[fid]->assets){
		this->makeDefault(fid);
		return this->nodes[fid]->assets / 2.0;
	}
	this->nodes[fid]->assets -= liquidate;
	// double sumL = 0;

	// for( auto & a: this->nodes[fid]->assets){
	// 	if(a.second + sumL > liquidate){
	// 		a.second = a.second - liquidate + sumL;
	// 		return 0;
	// 	}
	// 	else{
	// 		sumL += a.second;			
	// 		a.second = 0;
	// 	}
	// }
	// if(sumL<liquidate){
	// 	this->nodes[fid]->makeDefault();
	// 	// unwind(fid);
	// 	return sumL/2.0;
	// }
	return 0;
}

std::tuple<double,double,double> CreditNet::payAsset(int fid1, int fid2, double amt, string mode, int transSeqNum, string purpose, double drate, double haircut, double crate){
	this->updateNodeDegrees();
	if (mode == "SRC_DECIDE"){
		mode = this->nodes[fid1]->routePreference;
	}

	// this->print();
	if(verb){
		cout<<"purpose: "<<purpose<<endl;
		cout<<"amount: "<<amt<<endl;		
	}


	// fid1 = 2; fid2 = 0;
	// cout << "from " << fid1 << ", to " << fid2 << " " << request << " " << mode;

	CplexConverter converter;
	converter.constructCplex(this, this->nodes[fid1], this->nodes[fid2], amt, transSeqNum, verb);
	transactionCounter += 1;
	// converter.printInput();
	
	LpSolver lpSolver;

	this->nodes[fid1]->srcNum++;
	this->nodes[fid2]->destNum++;
	
	if (lpSolver.solveLpProblem(converter, mode, purpose, crate, deposit_rate, haircut)){
		
		double check = converter.copyBack();
		double debt_coc_rec = converter.getDebtCost();
		double asset_coc_rec = converter.getAssetCost();
		if (round(precision*check)/precision != round(amt*precision)/precision){
			// cout<< "unwinding amt check failed"<<endl;
			// cout<< "requested "<< amt<< " routed "<<check<<endl; 
			// converter.printResult();
			// cout<<"routed "<<check<<" out of "<<amt<<endl;
			// converter.printResult();
			return std::make_tuple(check,asset_coc_rec,debt_coc_rec);
		}
		// this->nodes[fid1]->transactionNum++;
		// this->nodes[fid1]->successSrc++;
		// this->nodes[fid2]->successDest++;
		// cout << " success " << endl;
		// converter.printResult();
		return std::make_tuple(0.0,asset_coc_rec,debt_coc_rec);
	}
	// cout << "infeasible payment" << endl;
	// if(amt < 0){
	// 	cout<<"negative requested amount: "<<amt<<endl;
	// }
	// converter.printResult();	
	return std::make_tuple(-1.0,0.0,0.0);
}

int CreditNet::pay(int fid1, int fid2, double amt, string mode, int transSeqNum){
	this->updateNodeDegrees();
	if (mode == "SRC_DECIDE"){
		mode = this->nodes[fid1]->routePreference;
	}

	// this->print();

	// fid1 = 2; fid2 = 0;
	// cout << "from " << fid1 << ", to " << fid2 << " " << request << " " << mode;

	CplexConverter converter;
	converter.constructCplex(this, this->nodes[fid1], this->nodes[fid2], amt, transSeqNum, verb);
	transactionCounter += 1;
	// converter.printInput();
	
	LpSolver lpSolver;

	this->nodes[fid1]->srcNum++;
	this->nodes[fid2]->destNum++;
	
	if (lpSolver.solveLpProblem(converter, mode, "DEBT",this->deposit_rate,this->deposit_rate,this->haircut)){
		double check = converter.copyBack();
		if (round(precision*check)/precision != round(amt*precision)/precision){
			// cout<< "unwinding amt check failed"<<endl;
			// cout<< "requested "<< amt<< " routed "<<check<<endl; 
			// converter.printResult();

			// return 1;
		}
		// cout<<"routed "<<check<<" out of "<<amt<<endl;
		// this->nodes[fid1]->transactionNum++;
		// this->nodes[fid1]->successSrc++;
		// this->nodes[fid2]->successDest++;
		// cout << " success " << endl;
		return 0;
	}
	cout << " fail " << endl;
	return 1;
}

// void CreditNet::unwind(int fid){


// 	// int numD = 0;
// 	// for (auto nodePair : nodes) {
// 	// if(nodes[fid]->getWealth() > 0){
// 		// this->nodes[fid]->print();	
// 		// cout<<"node "<<fid<<" before default collateral: "<<this->nodes[fid]->getWealth()<<endl;
// 		// cout<<"node "<<fid<<" before default wealth: "<<this->nodes[fid]->getScrip()<<endl;		
// 	// }
// 	grabCollateral(fid,this->nodes[fid]->sumAssets()* 0.5);
// 	double deposit_paid = payAsset(fid, marketId, this->nodes[fid]->deposits,"MAX_FLOW", transactionCounter++, "DEBT",this->deposit_rate,this->deposit_rate,this->haircut);
// 	double cashleft = this->nodes[fid]->getCash();
// 	std::vector<double> debt_pay(nodeNum -1, 0.0);
// 	double debtSum = 0;
// 	// Node* n = nodes[fid1];
// 	for (auto &eIn : this->nodes[fid]->edge_in){
// 		for (auto &cEdge : eIn.second->singleCreditEdges){
// 			for (auto &debtEdge : cEdge->debt_current){
// 				debt_pay[eIn.second->nodeFrom->nodeId] += debtEdge.second->capacity;
// 				debtSum += debtEdge.second->capacity;
// 			}
// 		}
// 	}

// 	for (int p = 0; p < nodeNum -1; p++){
// 		grabCollateral(debt_pay[p]/debtSum,p);
// 		payCollateral(fid,debt_pay[p]/debtSum);
// 		cReturns[fid] += debt_pay[p]/debtSum;
// 	}


// 	// if(deposit_paid == 0){
// 	// 	const int trxNodes = nodes.size() - 1;
// 	// 	std::vector<double> cr_pay(trxNodes, 0.0);
// 	// 	// Node* n = nodes[fid1];
// 	// 	for (auto &eIn : this->nodes[fid]->edge_in){
// 	// 		for (auto &cEdge : eIn.second->singleCreditEdges){
// 	// 			for (auto &debtEdge : cEdge->debt_current){
// 	// 				cr_pay[eIn.second->nodeFrom->nodeId] += debtEdge.second->capacity*cEdge->collateralRate;
// 	// 			}
// 	// 		}
// 	// 	}
// 	// 	std::vector<int> paid(trxNodes,0);
// 	// 	paid[fid] = 1;
// 	// 	bool payDone = false;
// 	// 	// while (not payDone){
// 	// 		bool payDone1 = true;
// 	// 		for (int i=0;i<cr_pay.size();++i){
// 	// 			if (paid[i] == 0 && i!= fid){
// 	// 				// cout<<"unwinding collateral "<<cr_pay[i]<<" to node "<<i<<endl;
// 	// 				double success = payAsset(fid, i, cr_pay[i], "MAX_FLOW", transactionCounter++,"DEBT",CR,deposit_rate,haircut);
// 	// 				if (success == 0){
// 	// 					// cout<<"collateral unwound"<<endl;
// 	// 					paid[i] = 1;
// 	// 					cReturns[i]+=cr_pay[i];
// 	// 					payDone1 = false;
// 	// 				}
// 	// 				else{
// 	// 					cReturns[i]+=success;
// 	// 					cr_pay[i] -= success;
// 	// 				}
// 	// 			}
// 	// 		}
// 	// 		// payDone = payDone1;
// 	// 	// }

// 	// 	std::vector<double> debt_pay(trxNodes, 0.0);
// 	// 	double debtSum = 0;
// 	// 	// Node* n = nodes[fid1];
// 	// 	for (auto &eIn : this->nodes[fid]->edge_in){
// 	// 		for (auto &cEdge : eIn.second->singleCreditEdges){
// 	// 			for (auto &debtEdge : cEdge->debt_current){
// 	// 				debt_pay[eIn.second->nodeFrom->nodeId] += debtEdge.second->capacity;
// 	// 				debtSum += debtEdge.second->capacity;
// 	// 			}
// 	// 		}
// 	// 	}

// 	// 	double remaining = 0;
// 	// 	for (auto &eOut : this->nodes[fid]->edge_out){
// 	// 		for (auto &cEdge : eOut.second->singleCreditEdges){
// 	// 			for (auto &debtEdge : cEdge->debt_current){
// 	// 				remaining += debtEdge.second->capacity;
// 	// 			}
// 	// 		}
// 	// 	}

// 	// 	std::vector<int> paidDebt(trxNodes,0);
// 	// 	paidDebt[fid] = 1;
// 	// 	bool debtDone = false;
// 	// 	while (not debtDone){
// 	// 		bool debtDone1 = true;
// 	// 		for (int i=0;i<trxNodes;++i){
// 	// 			if (paidDebt[i] == 0 && i!= fid){
// 	// 				// cout<<"unwinding debt "<<debt_pay[i]<<" to node "<<i<<endl;
// 	// 				double success = payAsset(fid, i, debt_pay[i] * remaining / (0.00000000001+debtSum), "MAX_FLOW", transactionCounter++,"DEBT",CR,deposit_rate,haircut);

// 	// 				if (success == 0){
// 	// 					paidDebt[i] = 1;
// 	// 					debtDone1 = false;
// 	// 					cReturns[i] += debt_pay[i] * remaining / (0.00000000001+debtSum);
// 	// 				}
// 	// 				else{
// 	// 					cReturns[i] += success;
// 	// 					debt_pay[i] -= success* debtSum/(0.00000000001+remaining);
// 	// 				}
// 	// 			}
// 	// 		}
// 	// 		debtDone = debtDone1;
// 	// 	}
// 	// 	// cout<<"after default collateral: "<<this->nodes[fid]->getWealth()<<endl;
// 	// 	// cout<<"after default wealth: "<<this->nodes[fid]->getScrip()<<endl;	
// 	// 	// this->nodes[fid]->print();
// 	// }
	
// }

void CreditNet::printPayoff(){
	for (auto it : nodes) {
		cout << it.second->getNodeId() << ": "
			<< "Transactions " << it.second->transactionNum << ", Current Balance "
			<< it.second->getCurrBalance() << endl;
	}
}

// for liquid stuff
int CreditNet::genInterBankTrans(double request, string mode, int transSeqNum){

	// Node* f1 = NULL;
	// Node* f2 = NULL;
	// int trxNodes = nodes.size()-2;

	// int fid1 = -1;
	// int fid2 = -1;

	// while (fid1 == fid2){
	// 	fid1 = credNetConstants.uniformIntDistribution(
	// 	credNetConstants.globalGenerator)%(trxNodes);		
	// 	// while (this->nodes[fid1]->defaulted){
	// 	// 	fid1 = credNetConstants.uniformIntDistribution(
	// 	// 	credNetConstants.globalGenerator)%(nodeNum - 1);
	// 	// }
	// 	fid2 = credNetConstants.uniformIntDistribution(
	// 	credNetConstants.globalGenerator)%(trxNodes);
	// 	// while (this->nodes[fid2]->defaulted){
	// 	// 	fid2 = credNetConstants.uniformIntDistribution(
	// 	// 	credNetConstants.globalGenerator)%(nodeNum - 1);
	// 	// }
	// }	


	
	// // int fid1 = rand()%nodeNum;
	// // int fid2 = rand()%nodeNum;
	// // while (fid1 == fid2){
	// // 	fid2 = rand()%nodeNum;
	// // }

	// this->updateNodeDegrees();
	// if (mode == "SRC_DECIDE"){
	// 	mode = this->nodes[fid1]->routePreference;
	// }

	// int collateralized = checkCollateral(fid1);
	// if (collateralized == 1){
	// 	// cout<<"Tried transaction to "<<fid2<<", source is leveraged"<<endl;
	// 	// this->nodes[fid1]->print();
	// 	return -1;
	// }

	// if (this->nodes[fid1]->defaulted){
	// 	// cout<<"Tried transaction, source is defaulted"<<endl;
	// 	return 1;
	// }
	// // this->print();

	// // fid1 = 2; fid2 = 0;
	// // cout << "from " << fid1 << ", to " << fid2 << " " << request << " " << mode;

	// CplexConverter converter;
	// converter.constructCplex(this, this->nodes[fid1], this->nodes[fid2], request, transSeqNum);
	// transactionCounter += 1;
	// // converter.printInput();
	
	// LpSolver lpSolver;

	// this->nodes[fid1]->srcNum++;
	// this->nodes[fid2]->destNum++;
	
	// if (lpSolver.solveLpProblem(converter, mode)){
	// 	double check = converter.copyBack();
	// 	if (round(precision*check)/precision != round(precision*request)/precision){
	// 		cout<< "transaction amt check failed from node "<<fid1<< " to node "<<fid2<<endl;
	// 		cout<< "requested " << request << " routed "<<check<<endl;
	// 		// converter.printResult();
			
	// 		// return 1;
	// 	}		
	// 	this->nodes[fid1]->transactionNum++;
	// 	this->nodes[fid1]->successSrc++;
	// 	this->nodes[fid2]->successDest++;
	// 	// cout << " success " << endl;
	// 	return 0;
	// }
	// cout << " fail " << endl;
	return 1;
}
int CreditNet::checkCollateral(int fid1){
	double col = this->nodes[fid1]->getCollateral(this->deposit_rate);
	// cout<<"got collateral"<<endl;
	double w = this->nodes[fid1]->getWealth(this->haircut);
	// this->nodes[fid1]->print();
	double buff = w - col;
	// cout<<"wealth: "<<w<<endl;
	// cout<<"collateral extra: "<<buff<<endl;
	// cout<<"defaulted? "<<this->nodes[fid1]->defaulted;
	if (buff >= 0){
		this->nodes[fid1]->makeLeveraged(false, this->haircut);		
		return 0;
	}
	// double target = buff/(-CR);
	// if(this->nodes[fid1]->getCash()>=target ){
	// 	payCollateral(fid1,target);
	// 	double sum = 0;
	// 	for (auto &it : this->nodes[fid1]->atomicEdge_out){
	// 		int fromId = it.second->nodeFrom->nodeId;			
	// 		if(it.second->isDebt && (sum + it.second->capacity)<target){
	// 			sum += it.second->capacity;
	// 			it.second->route(it.second->capacity, it.second->interest_rate, this, this->transactionCounter++);			
	// 			this->nodes[fid1]->creditPayOut[fromId].second = 0.0;
	// 		}
	// 		if(it.second->isDebt && (sum + it.second->capacity)>=target){
	// 			it.second->route(target - sum, it.second->interest_rate, this, this->transactionCounter++);
	// 			this->nodes[fid1]->makeLeveraged(false);
	// 			this->nodes[fid1]->creditPayOut[fromId].second = it.second->capacity * (1.0 + it.second->interest_rate);				
	// 			if(verb){
	// 				cout<<"unleveraged node "<<fid1<<endl;
	// 			}	
	// 			return 0;
	// 		}
	// 	}
	// 	if(sum < target){
	// 		this->nodes[fid1]->deposits -= (target - sum);
	// 		this->nodes[fid1]->makeLeveraged(false);
	// 		if(verb){
	// 			cout<<"unleveraged node "<<fid1<<endl;
	// 		}					
	// 		return 0;			
	// 	}
	// }
	// string mode = this->nodes[fid1]->routePreference;

	// CplexConverter converter;
	// converter.constructCplex(this, this->nodes[fid1], this->nodes[nodes.size()-1], -buff, transactionCounter);
	// transactionCounter += 1;
	// // converter.printInput();
	
	// LpSolver lpSolver;

	
	// if (lpSolver.solveLpProblem(converter, mode, "DEBT")){
	// 	// converter.printResult();
	// 	double check = converter.copyBack();
	// 	if (round(precision*check)/precision != round(precision*-buff)/precision){
	// 		cout<< "collateral amt check failed"<<endl;
	// 		cout<< "routed "<<check<< " needed "<< -buff <<endl;
	// 		// this->nodes[fid1]->makeLeveraged(true);
	// 		// return 1;
	// 	}		
	// 	// cout << "fcollateral success " << endl;		
	// 	grabCollateral(fid1,-buff);
	// 	this->nodes[fid1]->makeLeveraged(false);
	// 	// route buff from market to fid1
	// 	return 0;
	// }
	// else{
		// cout<<"collateral default "<<fid1<<endl;
	this->nodes[fid1]->makeLeveraged(true, this->haircut);
	return 1;
	// }
}

void CreditNet::grabCollateral(int fid, double amt){
	for (auto &it : this->nodes[nodes.size() - 1]->atomicEdge_in){
		if (it.second->nodeFrom->nodeId == fid){
			it.second->route(amt, 0, this, this->transactionCounter);			
		}
		transactionCounter += 1;
	}
}

void CreditNet::activateReserves(int fid){
	double amt = max(0.0,this->nodes[fid]->mReserved);
	
	grabCollateral(fid,amt);
	// for (auto &it : this->nodes[marketId]->atomicEdge_in){
	// 	if (it.second->nodeFrom->nodeId == fid){
	// 		it.second->route(amt, 0, this, this->transactionCounter);			
	// 	}
	// 	transactionCounter += 1;
	// }
	this->nodes[fid]->mReserved = 0.0;
	this->nodes[fid]->debt_reserve = 0.0;
}

void CreditNet::deactivateReserves(int fid, int mode){
	// min(this->nodes[fid]->reserveT,this->nodes[fid]->getCash());
	// cout<<amt<<endl;
	double amt = this->nodes[fid]->reserveT;
	if(mode == 1){
		//reserve for lending
		amt += this->nodes[fid]->debt_reserve;
	}
	amt = min(amt, this->nodes[fid]->getCash());

	double status = payCash(fid, amt);
	// for (auto &it : this->nodes[marketId]->atomicEdge_out){
	// 	if (it.second->nodeFrom->nodeId == fid){
	// 		it.second->route(amt, 0, this, this->transactionCounter);			
	// 	}
	// 	transactionCounter += 1;
	// }
	// cout<<"cash available "<<this->nodes[fid]->getCash()<<endl;
	this->nodes[fid]->mReserved = amt;
	// cout<<"deactivated"<<endl;
}


int CreditNet::payCollateral(int fid, double amt){
	for (auto &it : this->nodes[fid]->atomicEdge_in){
		if (it.second->nodeFrom->nodeId == marketId && it.second->isDebt){
			if(it.second->capacity < amt){
				// cout<<"cannot make cash payment"<<endl;
				return 1;
			}
			it.second->route(amt, 0, this, this->transactionCounter++);
			return 0;		
		}
	}
}

// static void deepCopyHelper(CreditNet* newGraph, CreditNet& oldGraph){

// 	newGraph->precision = oldGraph.precision;
// 	newGraph->transactionCounter = oldGraph.transactionCounter;
// 	newGraph->initVol = oldGraph.initVol;
// 	newGraph->returns = oldGraph.returns;
// 	newGraph->volatilities = oldGraph.volatilities;
// 	newGraph->wealths = oldGraph.wealths;
// 	newGraph->credits = oldGraph.credits;
// 	newGraph->credits_last = oldGraph.credits_last;
// 	newGraph->cReturns = oldGraph.cReturns;
// 	newGraph->rDefaults = oldGraph.rDefaults;
// 	newGraph->cDefaults = oldGraph.cDefaults;
// 	newGraph->dDefaults = oldGraph.dDefaults;
// 	newGraph->aDefaults = oldGraph.aDefaults;
// 	newGraph->ntrials = oldGraph.ntrials;
// 	newGraph->deposit_rate = oldGraph.deposit_rate;
// 	newGraph->haircut = oldGraph.haircut;
// 	newGraph->maturity = oldGraph.maturity;
// 	newGraph->initR = oldGraph.initR;
// 	newGraph->verb = oldGraph.verb;
// 	newGraph->fDefault = oldGraph.fDefault;

// 	newGraph->nodeNum = oldGraph.nodeNum;
// 	newGraph->atomicGlobalId = oldGraph.atomicGlobalId;
// 	newGraph->marketId= oldGraph.marketId;
// 	newGraph->expected_deposit= oldGraph.expected_deposit;
// 	newGraph->expected_asset_return= oldGraph.expected_asset_return;
// 	newGraph->asset_volatility= oldGraph.asset_volatility;
// 	newGraph->CR= oldGraph.CR;
// 	newGraph->deposit_shock= oldGraph.deposit_shock;
// 	newGraph->FFR= oldGraph.FFR;
// 	for (auto& it : oldGraph.nodes){
// 		newGraph->nodes[it.first] = new Node(it.first);
// 		newGraph->nodes[it.first]->transactionNum = it.second->transactionNum;
// 		newGraph->nodes[it.first]->routePreference = it.second->routePreference;
// 		newGraph->nodes[it.first]->srcNum = it.second->srcNum;
// 		newGraph->nodes[it.first]->destNum = it.second->destNum;
// 		newGraph->nodes[it.first]->successSrc = it.second->successSrc;
// 		newGraph->nodes[it.first]->successDest = it.second->successDest;
// 		newGraph->nodes[it.first]->transSeq = it.second->transSeq;
// 		newGraph->nodes[it.first]->degree = it.second->degree;
// 		newGraph->nodes[it.first]->defaulted = it.second->defaulted;
// 		newGraph->nodes[it.first]->leveraged = it.second->leveraged;
// 		newGraph->nodes[it.first]->isMarket = it.second->isMarket;
// 		newGraph->nodes[it.first]->theta = it.second->theta;
// 		newGraph->nodes[it.first]->deposits = it.second->deposits;
// 		newGraph->nodes[it.first]->lambda = it.second->lambda;
// 		newGraph->nodes[it.first]->w_assets = it.second->w_assets;
// 		newGraph->nodes[it.first]->assets = it.second->assets;
// 		newGraph->nodes[it.first]->creditReturn = it.second->creditReturn;
// 		newGraph->nodes[it.first]->creditVol = it.second->creditVol;
// 		newGraph->nodes[it.first]->credit_returns_in = it.second->credit_returns_in;
// 		newGraph->nodes[it.first]->credit_vol_in = it.second->credit_vol_in;
// 		newGraph->nodes[it.first]->folio_volume = it.second->folio_volume;
// 		newGraph->nodes[it.first]->creditPayIn = it.second->creditPayIn;
// 		newGraph->nodes[it.first]->creditPayOut = it.second->creditPayOut;
// 	}
// 	for (auto& oldNodePair : oldGraph.nodes){
// 		for (auto& oldEdgePair : oldNodePair.second->edge_in){
// 			int fromId = oldEdgePair.first;
// 			int toId = oldNodePair.first;
// 			Node* newFromNode = newGraph->nodes[fromId];
// 			Node* newToNode = newGraph->nodes[toId];

// 			Edge* e = new Edge(*(oldEdgePair.second),
// 				newFromNode, newToNode, newGraph->atomicEdges);
// 			newToNode->edge_in[fromId] = e;
// 			newFromNode->edge_out[toId] = e;
// 			newGraph->edges.push_back(e);
// 		}
// 	}

// }

// CreditNet::CreditNet(CreditNet &graphT){
// 	deepCopyHelper(this, graphT);
// }

// CreditNet& CreditNet::operator=(CreditNet &graphT){
// 	deepCopyHelper(this, graphT);
// 	return *this;
// }

// int CreditNet::payIR(int fid1){
// 	// // if (this->nodes[fid1]->defaulted){
// 	// 	// cout<<"node "<<fid1<< " is already defaulted"<<endl;
// 	// 	return 0;
// 	// }
// 	// // int numD = 0;
// 	// // for (auto nodePair : nodes) {
// 	// const int nodeSize = nodes.size();
// 	// std::vector<double> ir_pay(nodeSize, 0.0);
// 	// // Node* n = nodes[fid1];
// 	// for (auto &atoOut : this->nodes[fid1]->atomicEdge_out){
// 	// 	if ((atoOut.second->isDebt)&&(atoOut.second->capacity>0)){
// 	// 		// cout<<"debt detected: "<<atoOut.second->interest_rate * atoOut.second->capacity/100.0<<endl;
// 	// 		// atoOut.second->print();
// 	// 		int fid2 = atoOut.second->nodeFrom->getNodeId();
// 	// 		double rate = atoOut.second->interest_rate/100.0;
// 	// 		// cout<<"rate is "<<rate<<endl;
// 	// 		// cout<<"debt amount is "<<atoOut.second->capacity<<endl;
// 	// 		// cout<<atoOut.second->interest_rate * atoOut.second->capacity/100.0<<" amount interest paid to: "<<fid2<<" out of "<<nodeNum<<endl;
// 	// 		ir_pay[fid2] += rate * atoOut.second->capacity;
// 	// 		// ir_pay[fid2] += (int) atoOut.second->capacity/(10-atoOut.second->interest_rate);
// 	// 		// cout<<ir_pay[fid2]<<" will be paid to "<<fid2<<endl;
// 	// 	}
// 	// }
// 	// // for (int i=0; i< nodeSize; i++){
// 	// // 	checkCollateral(i);
// 	// // }	
// 	// // cout<< "added up "<<fid1<<endl;
// 	// for (int i=0; i< nodeSize; i++){
// 	// 	if (this->nodes[i]->defaulted || this->nodes[i]->isMarket){
// 	// 		// cout<<"node "<<i<< " is market or is already defaulted"<<endl;
// 	// 		continue;
// 	// 	}		
// 	// 	// cout<<ir_pay[i]<<" will be paid to "<<i<<endl;		
// 	// 	if(ir_pay[i]>0){
// 	// 		// cout<<"ir to pay by "<< fid1 << " to "<< i<< " is "<< ir_pay[i]<<" rounded to " << round(ir_pay[i]*precision)/precision<<endl;
// 	// 		// this->updateNodeDegrees();
// 	// 		// cout<<"updated degree"<<endl;
// 	// 		string mode = this->nodes[fid1]->routePreference;

// 	// 		// this->print();

// 	// 		// fid1 = 2; fid2 = 0;
// 	// 		// cout << "from " << fid1 << ", to " << fid2 << " " << request << " " << mode;

// 	// 		CplexConverter converter;
// 	// 		converter.constructCplex(this, this->nodes[fid1], this->nodes[i], round(ir_pay[i]*precision)/precision, transactionCounter);
// 	// 		transactionCounter += 1;
// 	// 		// converter.printInput();
			
// 	// 		LpSolver lpSolver;

			
// 	// 		if (lpSolver.solveLpProblem(converter, mode)){
// 	// 			// converter.printResult();
// 	// 			double check = converter.copyBack();
// 	// 			if (round(precision*check)/precision != round(precision*converter.request)/precision){
// 	// 				cout<< "interest amt check failed"<<endl;
// 	// 				cout<< "owed "<<ir_pay[i]<< " paid "<<check<<endl;
// 	// 				// cout<<"default"<<endl;
// 	// 				// this->nodes[fid1]->makeDefault();
// 	// 				// unwind(fid1);					
// 	// 				// return 1;
// 	// 			}				
// 	// 			// cout << " success " << endl;
// 	// 		}
// 	// 		// cout << " fail " << endl;
// 	// 		else{
// 	// 			// cout<<"default"<<endl;
// 	// 			// this->nodes[fid1]->print();	

// 	// 			// this->nodes[i]->print();	
// 	// 			// this->print();
// 	// 			this->nodes[fid1]->makeDefault();
// 	// 			unwind(fid1);				
// 	// 			// cout<<"deleted"<<endl;
// 	// 			// this->nodes[fid1]->print();
// 	// 			// cout<<"printed"<<endl;
// 	// 			return 1;
// 	// 		}
// 	// 	}
// 	// }
// 	template <typename T>
// 		vector<size_t> sort_indexes(const vector<T> &v) {

// 		  // initialize original index locations
// 		  	vector<size_t> idx(v.size());
// 		  	iota(idx.begin(), idx.end(), 0);

// 		  // sort indexes based on comparing values in v
// 			sort(idx.begin(), idx.end(),
// 		       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

// 		  return idx;
// 		}

// 	return 0;
// }
