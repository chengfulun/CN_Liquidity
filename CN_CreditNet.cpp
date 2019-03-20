// CN_CreditNet.C 
#define NON_DEBUG_MODE

#include "CN_Constants.h"
#include "CN_CreditNet.h"
#include <string>
#include <queue>
#include <vector>
#include <iostream>
#include <math.h>

using namespace std;

extern CredNetConstants credNetConstants;

CreditNet::CreditNet(int finNumT, double precision, int marketId, double initR,double initVol, double deposit_rate, double haircut, int maturity, bool verb) : Graph(finNumT, marketId){
	this->precision = precision;
	this->transactionCounter = 0;
	this->marketId = marketId;
	this->initVol = initVol;
	this->initR = initR;
	this->deposit_rate = deposit_rate;
	this->haircut = haircut;
	this->maturity = maturity;
	this->verb = verb;
	returns.assign(nodeNum, initR);
	wealths.assign(nodeNum, 0.0);
	credits.assign(nodeNum, 0.0);
	credits_last.assign(nodeNum, 0.0);
	rDefaults.assign(nodeNum, 0.0);
	cDefaults.assign(nodeNum, 0.0);	
	dDefaults.assign(nodeNum, 0.0);	
	aDefaults.assign(nodeNum, 0.0);	
	fDefaults.assign(nodeNum, 0.0);	
	liquid_loss.assign(nodeNum, 0.0);
	default_wealth.assign(nodeNum, 0.0);	
	deposit_losses.assign(nodeNum, 0.0);
	creditor_losses.assign(nodeNum, 0.0);	
	// credits_sq.assign(nodeNum, 0.0);

	cReturns.assign(nodeNum, 0.0);
	// cReturns_sq.assign(nodeNum, 0.0);
	// credits_cross.assign(nodeNum, 0.0);
	volatilities.assign(nodeNum, initVol);
	for (int i=0;i<nodeNum;i++){
		this->nodes[i]->defaulted = false;
	}
}

CreditNet::CreditNet() : Graph(){}

vector<double> CreditNet::credit_shares(int ix){
	vector<double> shares;
	vector<double> out;	
	double sum = 0;
	for (int i = 0; i<returns.size(); i++){
		if(not this->nodes[ix]->leveraged && i != ix && i != returns.size()-1 && not this->nodes[i]->defaulted){
			// cout<<"wealth "<<wealths[i]<<" returns "<<returns[i]<<" vols "<<volatilities[i]<<endl;

			// I think this is where sharpe ratio is calculated
			// double weight = max(0.0,wealths[i]*returns[i]/(sqrt(volatilities[i]))+0.00000000001);
			double pr_not_default = 1-this->nodes[i]->defaultProb();
			double var_not_default = pr_not_default * (1-pr_not_default);
			double weight = max(0.0, 0.00000000001+ pr_not_default * returns[i] /(
										var_not_default * (double)volatilities[i] 
										+ var_not_default * (double)returns[i] 
										+ pr_not_default * (double)returns[i]
									 )
								);

			// cout<<weight<<endl;
			shares.push_back(weight);
			sum += weight;
		}
		else{shares.push_back(0);}
	}
	// double sumOut = 0;
	// for (int j = 0; j<returns.size();j++){
	// 	double cand = shares[j]/sum;
	// 	if(cand < 0.1){
	// 		out.push_back(0);
	// 	}
	// 	else{
	// 		out.push_back(cand);
	// 		sumOut += cand;
	// 	}

	// }
	for (int j = 0; j<returns.size();j++){
		out.push_back(shares[j]/(0.0000000000001+sum));
	}
	return out;
}

vector<int> rix(int l){
	vector<int> v(l) ; // vector with 100 ints.
	iota (std::begin(v), std::end(v), 0); // Fill with 0, 1, ..., 99.
	random_shuffle ( v.begin(), v.end() );  // in place no extra array	
	return v;
}

int CreditNet::makeInvest(bool forced, bool verbose){

	int leverages = 0;
	updateReturns();
	// for(int i=0;i<returns.size() ; i++){
	// 	cout<<returns[i]<<endl;
	// // }

	// credits_last=credits;
	// update credits
	for(int cc = 0; cc<nodeNum; cc++){
		if (not this->nodes[cc]->defaulted){
			credits_last[cc] = credits[cc];			
			credits[cc] = 0.0;
		}
		if (this->nodes[cc]->defaulted && not this->nodes[cc]->been_defaulted){
			credits_last[cc] += credits[cc];
			cReturns[cc] += credits[cc];
		}
	}
	// fill(credits.begin(), credits.end(), 0.0);
	vector<int> randix  = rix(	nodeNum );
	for (int kk = 0; kk< nodeNum; kk++){
		int k = randix[kk];
		wealths[k] = this->nodes[k]->getWealth(haircut);
		leverages += checkCollateral(k);
		cDefaults[k] += checkCollateral(k);
		// cout<<"leveraged? "<<overleveraged<<endl;
		// if(this->nodes[j]->defaulted){
		// 	cout<<j<<" already defaulted"<<endl;
		// 	continue;
		// }
		if(forced && not this->nodes[k]->defaulted && this->fDefault > 0){
			this->nodes[k]->makeDefault();
			if(this->verb){
				cout<<"made "<<k<<" default "<<endl;				
			}
			this->fDefault -= 1;
			fDefaults[k] += 1;
		}
		if (k != marketId && not this->nodes[k]->defaulted){
			// cout<<j <<" not defaulted "<<this->nodes[j]->defaulted<<endl;
			// double w = wealths[j];
			int asset_maturity = maturity;
		
			//update debt payouts			

			// TODO calculate mu and sigma of credit
			double credit_mu = 0.0;
			double credit_sigma = 0.0;
			double c_weight = 0.0;
			double rsize = returns.size()-2.0; // number of (non-default ???) banks
			for (int ii = 0 ; ii < nodeNum; ii++){
				if(ii!=marketId && ii!=k){
					// this might not be right (the total credit I want) TODO
					double c_issue = credits_last[ii]; // credit issued to bank ii at the start of period
					
					double pr_not_default = 1-this->nodes[ii]->defaultProb();
					double var_not_default = pr_not_default * (1-pr_not_default);
					// for credit_mu calculation
					credit_mu += c_issue * pr_not_default * (double)returns[ii];
					// for credit_sigma calculation
					credit_sigma += c_issue * (
						var_not_default * (double)volatilities[ii] 
						+ var_not_default * (double)returns[ii] 
						+ pr_not_default * (double)returns[ii]
					);
					// add to total weight
					c_weight += c_issue; 
				}
			}
			if(c_weight == 0){
				credit_mu = initR;
				credit_sigma = initVol;
			}
			else{
				credit_mu = credit_mu / c_weight;
				credit_sigma = credit_sigma / c_weight;				
			}

			// calculate mu and sigma of credit
			// double c_average = 0.0; // credit return
			// double cv_average = 0.0; // credit volatility 
			// double wsum = 0.0;
			// double wsum2 = 0.0;			
			// double rsize = returns.size()-2.0; // number of (non-default ???) banks
			// for (int ii = 0 ; ii < nodeNum; ii++){
			// 	if(ii!=marketId && ii!=k){
			// 		// c_average += (double)returns[ii]/rsize;
			// 		// cv_average += (double)volatilities[ii]/rsize;
			// 		double c_issue = credits_last[ii]; // how much credit issued last period ???
			// 		// this->nodes[ii]->getCredit() + this->nodes[ii]->getDebt();
			// 		c_average += (double)returns[ii]*c_issue/rsize;
			// 		cv_average += (double)volatilities[ii]*c_issue*c_issue/rsize;
			// 		wsum += credits_last[ii];
			// 		wsum2 += credits_last[ii]*credits_last[ii];
			// 	}
			// }
			// if(wsum == 0){
			// 	c_average = initR;
			// 	cv_average = initVol;
			// }
			// else{
			// 	c_average = c_average / (wsum);
			// 	cv_average = cv_average / (  wsum2);				
			// }


			if(verb){
			//	cout<<"aggregate creturns "<<c_average<<" aggregate cvol "<<cv_average<<endl;
				cout<<"aggregate creturns "<<credit_mu<<" aggregate cvol "<<credit_sigma<<endl;
			}
			
			// accumulate( returns.begin(), returns.end(), 0.0)/returns.size();
			// if(c_average>0.02){
				// cout<<"credit returned "<<c_average<<endl;
			// }			
			 // = accumulate( volatilities.begin(), volatilities.end(), 0.0)/volatilities.size();

			// calculates lambda* and alpha* in section 3.4 of AI3 paper
			// this->nodes[k]->getLambda(this->expected_asset_return, this->asset_volatility*this->asset_volatility, c_average, cv_average, FFR, 1.0/(deposit_rate+0.00000000000000001));
			this->nodes[k]->getLambda(this->expected_asset_return, this->asset_volatility*this->asset_volatility, credit_mu, credit_sigma, FFR, 1.0/(deposit_rate+0.00000000000000001));
			if(verb){
				cout<<"node "<< k << " theta "<<this->nodes[k]->theta<<" portfolio size "<<this->nodes[k]->folio_volume<<" lambda "<<this->nodes[k]->lambda<<" assets "<<this->nodes[k]->w_assets<<" wealth "<<this->nodes[k]->getWealth(haircut)<<endl;
				cout<<"want assets "<<this->nodes[k]->getWealth(haircut)*this->nodes[k]->lambda*this->nodes[k]->w_assets<<endl;
			}

			// calculate actual asset goal
			double asset_target = this->nodes[k]->folio_volume*this->nodes[k]->w_assets;
			// if(overleveraged == 0){
			// calculate actual asset purchase amount (limit by paying ability)
			double fail = payAsset(k, this->marketId, asset_target, "MAX_FLOW", transactionCounter++, "ASSET",0.0,deposit_rate,haircut);
			if(verb){
				cout<<"paid assets "<<fail<<endl;
			}
			// this->nodes[i]->assets.push_back(std::make_pair(maturity,iAmount));

			// execute purchase	
			if(fail > 0){
				this->nodes[k]->buyAssets(fail,asset_maturity);
				// cout<<"asset shortage "<<asset_target - fail<<endl;
				// cout<<"got assets "<<failedil<<endl;
			}
			else{
				this->nodes[k]->buyAssets(asset_target,asset_maturity);
				// cout<<"got assets "<<this->nodes[j]->lambda*this->nodes[j]->w_assets*w<<endl;				
			}
		}	
	}

	for (int jjj = 0; jjj < nodeNum -1; jjj++){
		// this->nodes[jjj]->print();
		if(not this->nodes[jjj]->defaulted){
			cReturns[jjj] = this->nodes[jjj]->getCredit(); // incoming credit of nodes
		}
	}
	if(verbose){
		this->resultsOut();
	}
	for (int j = 0; j < nodeNum - 1; j++){

		// for non-default banks
		// clear all the debt it lend out to non-default banks
		// route something ???
		if(not this->nodes[j]->defaulted){
			for(auto& it : this->nodes[j]->atomicEdge_out){
				int fromId = it.second->nodeFrom->nodeId;
				if(it.second->isDebt && it.second->capacity>0 && not this->nodes[fromId]->defaulted){
					unordered_map<int,pair<double,double>>::const_iterator got = this->nodes[j]->creditPayOut.find (fromId);
					int tempix = it.second->singleCreditIndex;
					double tempCR = it.second->originEdge->singleCreditEdges[tempix]->collateralRate;
					if ( got == this->nodes[j]->creditPayOut.end() ){
						this->nodes[j]->creditPayOut[fromId] = make_pair(tempCR,it.second->capacity*(1.0+it.second->interest_rate));
					}
					else{
						// cout<<"before "<<this->nodes[j]->creditPayOut[fromId].second<<endl;
						this->nodes[j]->creditPayOut[fromId].second += it.second->capacity*(1.0+it.second->interest_rate);
						// cout<<"after "<<this->nodes[j]->creditPayOut[fromId].second<<endl;
					}
					// cout<<"cpayout "<<this->nodes[j]->creditPayOut[fromId].second<<" from "<<j<<" to "<<fromId<<endl;

					it.second->route(it.second->capacity,it.second->interest_rate,this,this->transactionCounter++);
				}

			}			
		}

		// for banks just defaulted
		// clear its debt ???
		// route that somewhere ???
		if( this->nodes[j]->defaulted && not this->nodes[j]->been_defaulted){
			for(auto& it : this->nodes[j]->atomicEdge_out){
				int fromId = it.second->nodeFrom->nodeId;

				if(it.second->isDebt && it.second->capacity>0){
					cReturns[j] -= it.second->capacity;
						// cout<<"after "<<this->nodes[j]->creditPayOut[fromId].second<<endl;
				}
					// cout<<"cpayout "<<this->nodes[j]->creditPayOut[fromId].second<<" from "<<j<<" to "<<fromId<<endl;

				it.second->route(it.second->capacity,it.second->interest_rate,this,this->transactionCounter++);
			}
			this->nodes[j]->been_defaulted = true;

		}					
	}	
	

	
	

	for (int j = 0; j < nodeNum -1; j++){
		// cout<<"leveraged? "<<overleveraged<<endl;
		// if(this->nodes[j]->defaulted){
		// 	cout<<j<<" already defaulted"<<endl;
		// 	continue;
		// }

		if (j != marketId && not this->nodes[j]->defaulted ){
			// cout<<j <<" not defaulted "<<this->nodes[j]->defaulted<<endl;&&not this->nodes[j]->leveraged

			// double w = wealths[j];
			// int overleveraged = checkCollateral(j);

			// }



			vector<double> cShares = credit_shares(j);
			// double ir = (credNetConstants.uniformIntDistribution(
			// 	credNetConstants.globalGenerator) % numIR);
			// cout<<"random: "<<ir<<endl;	
				// this->addUnitEdge(nodes.find(i)->second, nodes.find(j)->second, ir, rand()%2);
			for (int i = 0; i< nodeNum; i++){
			// cout<<cShares[i]<<endl;
				if(i!= j && i!= marketId){

					double cap = cShares[i] * this->nodes[j]->folio_volume*(1.0-this->nodes[j]->w_assets);
					// if(cap<0){cout<<cap<< " is credit"<<endl;}
					this->modifyCredit(j, i, max(0.0,cap));
					// this->nodes[i]->print();
					// credits_sq[i] += cap*cap;

					//cout<<credits[i]<<" "; //output credit line
					//cout.flush(); // output credit line
				}else{
					//cout<<"X "; // output credit line
				}
			}				
			//cout<<endl; // output credit line
			


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

	for (int i = 0; i< nodeNum - 1; i++){
		// cout<<"c is "<<credits[i]<<endl;
		// cout<<"cl is "<<credits_last[i]<<endl;
		// cout<<"sumcredit is "<<this->nodes[i]->getCredit()<<endl;
		if(not this->nodes[i]->defaulted){
			credits[i] += this->nodes[i]->getCredit();			
		}

		// this->nodes[i]->print();
	}
	return leverages;	
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
		cout<<ntrials<<" "<<this->nodes[i]->theta<< " "<<rDefaults[i]
		<<" "<<aDefaults[i]<<" "<<dDefaults[i]<<" "<<cDefaults[i]<<" "<<this->nodes[i]->getWealth(haircut)
		<<" "<<this->nodes[i]->getDebt()<<" "<<
		this->nodes[i]->getCredit()<<" "<<
		this->nodes[i]->getScrip() - this->nodes[i]->getCash()<<" "<<
		this->nodes[i]->folio_volume*(1-this->nodes[i]->w_assets)<<" "<<this->nodes[i]->deposits<<" "<<this->nodes[i]->getCash()<< " "<<this->nodes[i]->sumAssets()<<endl;
	}
}

void CreditNet::resultsOut_1(){
	// cout<<"printing"<<endl;
	for (int i = 0; i< nodeNum - 1; i++){
		cout<<ntrials<<" "<<this->nodes[i]->theta<< " "<<rDefaults[i]
		<<" "<<aDefaults[i]<<" "<<dDefaults[i]<<" "<<cDefaults[i]<<" "<<this->nodes[i]->getWealth(haircut)
		<<" "<<this->nodes[i]->getDebt()<<" "<<
		this->nodes[i]->getCredit()<<" "<<
		this->nodes[i]->getScrip() - this->nodes[i]->getCash()<<" "<<
		this->nodes[i]->folio_volume*(1-this->nodes[i]->w_assets)<<" "<<alternate<<" "<<fDefaults[i]<< " "<< this->nodes[i]->deposits<<" "<<this->nodes[i]->getCash()<< " "<<this->nodes[i]->sumAssets()<<
		" "<<liquid_loss[i]<<" "<<default_wealth[i]<< " "<<deposit_losses[i]<<" "<<creditor_losses[i]<<endl;
	}
}

int CreditNet::shockPay(double alpha, bool crisis){
	// fill(cReturns.begin(), cReturns.end(), 0.0);
	// cReturns[fid] = 0.0;	
	// if(this->nodes[fid]->defaulted){
	// 	cout<<fid<<" already defaulted"<<endl;
	// 	return 0;
	// }
	// cout<<fid <<" not defaulted shock "<<this->nodes[fid]->defaulted<<endl;
	int defaultC = 0;
	for (int fid = 0;fid<nodeNum - 1;fid++){
		if(fid!=marketId && not this->nodes[fid]->defaulted){
			int aSize = this->nodes[fid]->assets.size();
			for (int i = 0; i<this->nodes[fid]->assets.size();i++){

				// cout<<"maturity is "<<this->nodes[fid]->assets[i].first<<" amount is "<<this->nodes[fid]->assets[i].second<<endl;			
				// cout<<"index "<<i<<" size "<<aSize<<endl;
				this->nodes[fid]->assets[i].first = this->nodes[fid]->assets[i].first - 1;
				// cout<<"maturity is "<<this->nodes[fid]->assets[i].first<<" amount is "<<this->nodes[fid]->assets[i].second<<endl;
				if ( this->nodes[fid]->assets[i].first == 0){
					
					// cout<<"alpha is "<<alpha<<endl;
					grabCollateral(this->nodes[fid]->nodeId, this->nodes[fid]->assets[i].second*(alpha));
					// cout<<"asset paid off "<<this->nodes[fid]->assets[i].second*(alpha)<<endl;
				}
			}
			int jj = 0;
			while(jj != this->nodes[fid]->assets.size()){
				if (this->nodes[fid]->assets[jj].first == 0){
					this->nodes[fid]->assets.erase(this->nodes[fid]->assets.begin() + jj);
				}
				else{jj++;}
			}
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
	vector<int>randix  = rix(	nodeNum );

	for (int fid1 = 0;fid1<nodeNum;fid1++){
		int fid = randix[fid1];
		if(fid!=marketId && not this->nodes[fid]->defaulted){
			double dShock = ( (credNetConstants.uniformDoubleDistribution(
								credNetConstants.globalGenerator)) - 0.5) * 2.0* deposit_shock * this->nodes[fid]->getWealth(1.0);
			// cout<<"dShock is "<<dShock<<endl;
			if(not crisis || (crisis && this->cCount <= 0)){
				// double dShock = (credNetConstants.normalDistribution(
				// 								credNetConstants.globalGenerator)) * deposit_shock * this->nodes[fid]->deposits;
				// double dShock = ( (credNetConstants.uniformDoubleDistribution(
				// 					credNetConstants.globalGenerator)) - 0.5) * 2.0* deposit_shock * this->nodes[fid]->getWealth(1.0);
				
				this->nodes[fid]->deposits += dShock;
				if(this->nodes[fid]->deposits<0){
					dShock = dShock - this->nodes[fid]->deposits;
					this->nodes[fid]->deposits = 0.0;
				}
				if(verb){
					cout<<"node "<<fid<<" deposit shock "<<dShock<<" cash "<<this->nodes[fid]->getCash()<< " credit "<<this->nodes[fid]->getCredit()<<" assets "<< this->nodes[fid]->sumAssets()<<" leveraged "<<this->nodes[fid]->leveraged<<endl;
				}
				if(dShock >= 0){
					grabCollateral(fid,dShock);
				}
				else{
					double failedDeposit = payAsset(fid, this->marketId, -dShock, "MAX_FLOW", transactionCounter++, "DEBT",deposit_rate, deposit_rate, haircut);
					// cout<<"pay deposit"<<endl;
					if (failedDeposit > 0){
						// cout<<"was defaulted "<< this->nodes[fid]->defaulted<<endl;
							double status = makeLiquidate(fid,-dShock - failedDeposit);
							liquid_loss[fid] += status;

							if(status >0){
								if(verb){
									cout<<fid << " couldn't pay " << -dShock - failedDeposit<< " is defaulted "<<this->nodes[fid]->defaulted<<endl;								
								}
								// this->nodes[fid]->makeDefault();
								defaultC+=1;
								dDefaults[fid] +=1;
								double dlost = this->nodes[fid]->deposits - this->nodes[fid]->getCash() +(-dShock - failedDeposit) -status;
								default_wealth[fid] += max(0.0,this->nodes[fid]->getWealth(1.0));
								deposit_losses[fid] += max(0.0,dlost);
								creditor_losses[fid] += min(this->nodes[fid]->getDebt()+dlost,this->nodes[fid]->getDebt());
						}

						
					}
				}				
			}
			if (crisis && this->cCount>0){
				if(verb){
					cout<<"crisis at "<<fid<<endl;
				}
				if(verb){
					cout<<"node "<<fid<<" deposit shock "<<this->nodes[fid]->deposits<<" cash "<<this->nodes[fid]->getCash()<< " credit "<<this->nodes[fid]->getCredit()<<" assets "<< this->nodes[fid]->sumAssets()<<endl;
				}				
				this->cCount += -1;
				fDefaults[fid] +=1;										
				double failedDeposit = payAsset(fid, this->marketId, this->nodes[fid]->deposits, "MAX_FLOW", transactionCounter++, "DEBT",deposit_rate, deposit_rate, haircut);
				if (failedDeposit > 0){
					// cout<<"was defaulted "<< this->nodes[fid]->defaulted<<endl;
					double status = makeLiquidate(fid,this->nodes[fid]->deposits - failedDeposit);
					liquid_loss[fid] += status;
					if(status >0){
						if(verb){
							cout<<fid << " couldn't pay " << this->nodes[fid]->deposits - failedDeposit<< " is defaulted "<<this->nodes[fid]->defaulted<<endl;								
						}
						// this->nodes[fid]->makeDefault();
						defaultC+=1;
						dDefaults[fid] +=1;	
						double dlost = this->nodes[fid]->deposits - this->nodes[fid]->getCash() +(this->nodes[fid]->deposits - failedDeposit) -status;
						default_wealth[fid] += max(0.0,this->nodes[fid]->getWealth(1.0));
						deposit_losses[fid] += max(0.0,dlost);
						creditor_losses[fid] += min(this->nodes[fid]->getDebt()+dlost,this->nodes[fid]->getDebt());											

					}
			
				}
				this->nodes[fid]->deposits = 0;				
			}
		}
	}
	randix  = rix(	nodeNum );

	for (int f = 0;f<nodeNum ;f++){
		int ff = randix[f];
		if( ff != marketId && not this->nodes[ff]->defaulted){
			for (auto& cPay : this->nodes[ff]->creditPayOut){
				// cout<<"pay credit"<<endl;
				double current_crate = cPay.second.first ;
				// cout<<"c to pay "<<cPay.second.second<<" from "<<ff<<" to "<<cPay.first<<endl;
				double paymentFail = payAsset(ff, cPay.first, cPay.second.second,"MAX_FLOW", transactionCounter++, "DEBT",current_crate, deposit_rate, haircut);
				if(paymentFail > 0){
					double status = makeLiquidate(ff,cPay.second.second - paymentFail);
					if(status >0){
						// this->nodes[fid]->makeDefault();
						if(verb){
							cout<<"node "<<ff<<" couldn't pay credit "<<cPay.second.second - paymentFail<<" defaulted"<<endl;	
						}
						cPay.second.second = 0;
						defaultC+=1;
						rDefaults[ff] +=1;						

						grabCollateral(cPay.first,status);
						cReturns[ff] += status+paymentFail;											
					}
					else{
						grabCollateral(cPay.first,cPay.second.second - paymentFail);
						cReturns[ff] += cPay.second.second;
					}
				}
				else{
					cReturns[ff] += cPay.second.second;					
				}
				// cReturns_sq[fid] += cPay.second.second*cPay.second.second;
				// credits_cross[fid] += cPay.second.second * this->cap;
				cPay.second.second = 0;
				// cout<<"cpaid total "<<cReturns[ff]<<endl;

			}
			// if(this->nodes[ff]->defaulted){
			// 	this->nodes[ff]->been_defaulted += 1;
			// 	credits_last[ff]+= credits[ff];
			// 	cReturns[ff] += credits[ff] - ;
			// 	// returns[ff] = (cReturns[ff] / credits_last[ff]) - 1.0;
			// }
		}		
	}
	for(int dd = 0; dd < nodeNum - 1; dd++){
		if(not this->nodes[dd]->defaulted && this->nodes[dd]->getWealth(haircut) < 0.0){
			this->nodes[dd]->makeDefault();
			defaultC+=1;
			aDefaults[dd] +=1;
			double dlost = this->nodes[dd]->deposits - this->nodes[dd]->getCash() - (0.5*this->nodes[dd]->sumAssets());
			liquid_loss[dd] += (0.5*this->nodes[dd]->sumAssets());
			deposit_losses[dd] += max(0.0,dlost);
			creditor_losses[dd] += min(this->nodes[dd]->getDebt()+dlost,this->nodes[dd]->getDebt());
			// unwind(dd);
			if(verb){
				cout<<dd<<" negative wealth couldn't ... default "<<this->nodes[dd]->getWealth()<<endl;
			}

		}	
	}
	if(verb){
		cout<<"loop done"<<endl;
	}

	return defaultC;		
}

/*
 * Update variables: returns, volatilities 
 */
void CreditNet::updateReturns(){
	if(ntrials <3){
		for(int i = 0; i<nodeNum;i++){
			if(not this->nodes[i]->defaulted){ 
				// if not default, in first 3 periods, keep returns & volatilities at initial value
				returns[i] = this->initR;
				volatilities[i] = initVol;				
			}
		}
	}
	else{
		for (int i = 0; i< nodeNum -1; i++){
			double delta;
			double delta2;
			// if(this->nodes[i]->defaulted){
			// 	returns[i] = 0.0;
			// 	volatilities[i] = 0.0;
			// }
			// else{
				if(credits_last[i]<0.00001){
					delta = 0.0;
				}
				// else if(aDefaults[i] == 1){
				// 	delta = -returns[i];
				// }
				else{
					delta = (cReturns[i] - credits_last[i])/(credits_last[i] + 0.00000000000001) - returns[i];
				}
				// if(cReturns[i] > 0){
					// cout<<"credit returns "<<cReturns[i]<<" credit expenditure "<<credits_last[i]<<endl;
				// }
				if(cReturns[i]>credits_last[i]*1.1 && verb){
					cout<<"credits violated "<<i<<endl;
				}
				// cout<<"leveraged "<<this->nodes[i]->leveraged<<" defaulted "<<this->nodes[i]->defaulted<<endl;
				returns[i] = returns[i] + delta/(double) ntrials;
				// returns[i] = cReturns[i]/credits[i];
				if(verb){
					cout<<"node "<<i<<" returns delta "<<delta<<" returns "<<cReturns[i]<<" credits "<<credits_last[i]<<" defaulted "<<this->nodes[i]->defaulted<<endl;
				}
				// if(this->nodes[i]->defaulted){
				// 	delta2 = -
				// }
				delta2 = (cReturns[i] - credits_last[i])/(credits_last[i] + 0.00000000000001) - returns[i];
				// cout<<"vols delta "<<delta2<<endl;

				volatilities[i] = (volatilities[i]*(double) (ntrials-2) + (delta*delta2))/(double) (ntrials - 1);

			// }
		}
	}
	ntrials += 1;
}

double CreditNet::makeLiquidate(int fid, double amt){
	if(verb){
		cout<<"liquidating "<<amt<<" for node "<<fid<<endl;
	}
	double liquidate = amt * 2.0;
	double sumL = 0;
	for( auto & a: this->nodes[fid]->assets){
		if(a.second + sumL > liquidate){
			a.second = a.second - liquidate + sumL;
			return 0;
		}
		else{
			sumL += a.second;			
			a.second = 0;
		}
	}
	if(sumL<liquidate){
		this->nodes[fid]->makeDefault();
		// unwind(fid);
		return sumL/2.0;
	}
	return 0;
}

double CreditNet::payAsset(int fid1, int fid2, double amt, string mode, int transSeqNum, string purpose, double crate, double drate, double haircut){
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
	
	if (lpSolver.solveLpProblem(converter, mode, purpose, crate, deposit_rate, haircut)){
		double check = converter.copyBack();
		if (round(precision*check)/precision != round(amt*precision)/precision){
			// cout<< "unwinding amt check failed"<<endl;
			// cout<< "requested "<< amt<< " routed "<<check<<endl; 
			// converter.printResult();
			// cout<<"routed "<<check<<" out of "<<amt<<endl;

			return check;
		}
		// this->nodes[fid1]->transactionNum++;
		// this->nodes[fid1]->successSrc++;
		// this->nodes[fid2]->successDest++;
		// cout << " success " << endl;
		return 0;
	}
	// cout << " fail " << endl;
	return -1;
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

void CreditNet::unwind(int fid){


	// int numD = 0;
	// for (auto nodePair : nodes) {
	// if(nodes[fid]->getWealth() > 0){
		// this->nodes[fid]->print();	
		// cout<<"node "<<fid<<" before default collateral: "<<this->nodes[fid]->getWealth()<<endl;
		// cout<<"node "<<fid<<" before default wealth: "<<this->nodes[fid]->getScrip()<<endl;		
	// }
	grabCollateral(fid,this->nodes[fid]->sumAssets()* 0.5);
	double deposit_paid = payAsset(fid, marketId, this->nodes[fid]->deposits,"MAX_FLOW", transactionCounter++, "DEBT",this->deposit_rate,this->deposit_rate,this->haircut);
	double cashleft = this->nodes[fid]->getCash();
	std::vector<double> debt_pay(nodeNum -1, 0.0);
	double debtSum = 0;
	// Node* n = nodes[fid1];
	for (auto &eIn : this->nodes[fid]->edge_in){
		for (auto &cEdge : eIn.second->singleCreditEdges){
			for (auto &debtEdge : cEdge->debt_current){
				debt_pay[eIn.second->nodeFrom->nodeId] += debtEdge.second->capacity;
				debtSum += debtEdge.second->capacity;
			}
		}
	}

	for (int p = 0; p < nodeNum -1; p++){
		grabCollateral(debt_pay[p]/debtSum,p);
		payCollateral(fid,debt_pay[p]/debtSum);
		cReturns[fid] += debt_pay[p]/debtSum;
	}


	// if(deposit_paid == 0){
	// 	const int trxNodes = nodes.size() - 1;
	// 	std::vector<double> cr_pay(trxNodes, 0.0);
	// 	// Node* n = nodes[fid1];
	// 	for (auto &eIn : this->nodes[fid]->edge_in){
	// 		for (auto &cEdge : eIn.second->singleCreditEdges){
	// 			for (auto &debtEdge : cEdge->debt_current){
	// 				cr_pay[eIn.second->nodeFrom->nodeId] += debtEdge.second->capacity*cEdge->collateralRate;
	// 			}
	// 		}
	// 	}
	// 	std::vector<int> paid(trxNodes,0);
	// 	paid[fid] = 1;
	// 	bool payDone = false;
	// 	// while (not payDone){
	// 		bool payDone1 = true;
	// 		for (int i=0;i<cr_pay.size();++i){
	// 			if (paid[i] == 0 && i!= fid){
	// 				// cout<<"unwinding collateral "<<cr_pay[i]<<" to node "<<i<<endl;
	// 				double success = payAsset(fid, i, cr_pay[i], "MAX_FLOW", transactionCounter++,"DEBT",CR,deposit_rate,haircut);
	// 				if (success == 0){
	// 					// cout<<"collateral unwound"<<endl;
	// 					paid[i] = 1;
	// 					cReturns[i]+=cr_pay[i];
	// 					payDone1 = false;
	// 				}
	// 				else{
	// 					cReturns[i]+=success;
	// 					cr_pay[i] -= success;
	// 				}
	// 			}
	// 		}
	// 		// payDone = payDone1;
	// 	// }

	// 	std::vector<double> debt_pay(trxNodes, 0.0);
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

	// 	double remaining = 0;
	// 	for (auto &eOut : this->nodes[fid]->edge_out){
	// 		for (auto &cEdge : eOut.second->singleCreditEdges){
	// 			for (auto &debtEdge : cEdge->debt_current){
	// 				remaining += debtEdge.second->capacity;
	// 			}
	// 		}
	// 	}

	// 	std::vector<int> paidDebt(trxNodes,0);
	// 	paidDebt[fid] = 1;
	// 	bool debtDone = false;
	// 	while (not debtDone){
	// 		bool debtDone1 = true;
	// 		for (int i=0;i<trxNodes;++i){
	// 			if (paidDebt[i] == 0 && i!= fid){
	// 				// cout<<"unwinding debt "<<debt_pay[i]<<" to node "<<i<<endl;
	// 				double success = payAsset(fid, i, debt_pay[i] * remaining / (0.00000000001+debtSum), "MAX_FLOW", transactionCounter++,"DEBT",CR,deposit_rate,haircut);

	// 				if (success == 0){
	// 					paidDebt[i] = 1;
	// 					debtDone1 = false;
	// 					cReturns[i] += debt_pay[i] * remaining / (0.00000000001+debtSum);
	// 				}
	// 				else{
	// 					cReturns[i] += success;
	// 					debt_pay[i] -= success* debtSum/(0.00000000001+remaining);
	// 				}
	// 			}
	// 		}
	// 		debtDone = debtDone1;
	// 	}
	// 	// cout<<"after default collateral: "<<this->nodes[fid]->getWealth()<<endl;
	// 	// cout<<"after default wealth: "<<this->nodes[fid]->getScrip()<<endl;	
	// 	// this->nodes[fid]->print();
	// }
	
}

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

/**
 * check if over-leveraged and update status related to leverage
 * w: wealth
 * col: total sum collateral requirements
 */
int CreditNet::checkCollateral(int fid1){
	double col = this->nodes[fid1]->getCollateral(this->deposit_rate);
	// cout<<"got collateral"<<endl;
	double w = this->nodes[fid1]->getWealth(this->haircut);
	double buff = w - col;
	// cout<<"wealth: "<<w<<endl;
	// cout<<"collateral short: "<<buff<<endl;
	if (buff >= 0){
		this->nodes[fid1]->makeLeveraged(false);		
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
	// 	// cout << "collateral success " << endl;		
	// 	grabCollateral(fid1,-buff);
	// 	this->nodes[fid1]->makeLeveraged(false);
	// 	// route buff from market to fid1
	// 	return 0;
	// }
	// else{
		// cout<<"collateral default "<<fid1<<endl;
		this->nodes[fid1]->makeLeveraged(true);
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

static void deepCopyHelper(CreditNet* newGraph, CreditNet& oldGraph){

	newGraph->precision = oldGraph.precision;
	newGraph->transactionCounter = oldGraph.transactionCounter;
	newGraph->initVol = oldGraph.initVol;
	newGraph->returns = oldGraph.returns;
	newGraph->volatilities = oldGraph.volatilities;
	newGraph->wealths = oldGraph.wealths;
	newGraph->credits = oldGraph.credits;
	newGraph->credits_last = oldGraph.credits_last;
	newGraph->cReturns = oldGraph.cReturns;
	newGraph->rDefaults = oldGraph.rDefaults;
	newGraph->cDefaults = oldGraph.cDefaults;
	newGraph->dDefaults = oldGraph.dDefaults;
	newGraph->aDefaults = oldGraph.aDefaults;
	newGraph->ntrials = oldGraph.ntrials;
	newGraph->deposit_rate = oldGraph.deposit_rate;
	newGraph->haircut = oldGraph.haircut;
	newGraph->maturity = oldGraph.maturity;
	newGraph->initR = oldGraph.initR;
	newGraph->verb = oldGraph.verb;
	newGraph->fDefault = oldGraph.fDefault;

	newGraph->nodeNum = oldGraph.nodeNum;
	newGraph->atomicGlobalId = oldGraph.atomicGlobalId;
	newGraph->marketId= oldGraph.marketId;
	newGraph->expected_deposit= oldGraph.expected_deposit;
	newGraph->expected_asset_return= oldGraph.expected_asset_return;
	newGraph->asset_volatility= oldGraph.asset_volatility;
	newGraph->CR= oldGraph.CR;
	newGraph->deposit_shock= oldGraph.deposit_shock;
	newGraph->FFR= oldGraph.FFR;
	for (auto& it : oldGraph.nodes){
		newGraph->nodes[it.first] = new Node(it.first);
		newGraph->nodes[it.first]->transactionNum = it.second->transactionNum;
		newGraph->nodes[it.first]->routePreference = it.second->routePreference;
		newGraph->nodes[it.first]->srcNum = it.second->srcNum;
		newGraph->nodes[it.first]->destNum = it.second->destNum;
		newGraph->nodes[it.first]->successSrc = it.second->successSrc;
		newGraph->nodes[it.first]->successDest = it.second->successDest;
		newGraph->nodes[it.first]->transSeq = it.second->transSeq;
		newGraph->nodes[it.first]->degree = it.second->degree;
		newGraph->nodes[it.first]->defaulted = it.second->defaulted;
		newGraph->nodes[it.first]->leveraged = it.second->leveraged;
		newGraph->nodes[it.first]->isMarket = it.second->isMarket;
		newGraph->nodes[it.first]->theta = it.second->theta;
		newGraph->nodes[it.first]->deposits = it.second->deposits;
		newGraph->nodes[it.first]->lambda = it.second->lambda;
		newGraph->nodes[it.first]->w_assets = it.second->w_assets;
		newGraph->nodes[it.first]->assets = it.second->assets;
		newGraph->nodes[it.first]->creditReturn = it.second->creditReturn;
		newGraph->nodes[it.first]->creditVol = it.second->creditVol;
		newGraph->nodes[it.first]->credit_returns_in = it.second->credit_returns_in;
		newGraph->nodes[it.first]->credit_vol_in = it.second->credit_vol_in;
		newGraph->nodes[it.first]->folio_volume = it.second->folio_volume;
		newGraph->nodes[it.first]->creditPayIn = it.second->creditPayIn;
		newGraph->nodes[it.first]->creditPayOut = it.second->creditPayOut;
	}
	for (auto& oldNodePair : oldGraph.nodes){
		for (auto& oldEdgePair : oldNodePair.second->edge_in){
			int fromId = oldEdgePair.first;
			int toId = oldNodePair.first;
			Node* newFromNode = newGraph->nodes[fromId];
			Node* newToNode = newGraph->nodes[toId];

			Edge* e = new Edge(*(oldEdgePair.second),
				newFromNode, newToNode, newGraph->atomicEdges);
			newToNode->edge_in[fromId] = e;
			newFromNode->edge_out[toId] = e;
			newGraph->edges.push_back(e);
		}
	}

}

CreditNet::CreditNet(CreditNet &graphT){
	deepCopyHelper(this, graphT);
}

CreditNet& CreditNet::operator=(CreditNet &graphT){
	deepCopyHelper(this, graphT);
	return *this;
}

int CreditNet::payIR(int fid1){
	// // if (this->nodes[fid1]->defaulted){
	// 	// cout<<"node "<<fid1<< " is already defaulted"<<endl;
	// 	return 0;
	// }
	// // int numD = 0;
	// // for (auto nodePair : nodes) {
	// const int nodeSize = nodes.size();
	// std::vector<double> ir_pay(nodeSize, 0.0);
	// // Node* n = nodes[fid1];
	// for (auto &atoOut : this->nodes[fid1]->atomicEdge_out){
	// 	if ((atoOut.second->isDebt)&&(atoOut.second->capacity>0)){
	// 		// cout<<"debt detected: "<<atoOut.second->interest_rate * atoOut.second->capacity/100.0<<endl;
	// 		// atoOut.second->print();
	// 		int fid2 = atoOut.second->nodeFrom->getNodeId();
	// 		double rate = atoOut.second->interest_rate/100.0;
	// 		// cout<<"rate is "<<rate<<endl;
	// 		// cout<<"debt amount is "<<atoOut.second->capacity<<endl;
	// 		// cout<<atoOut.second->interest_rate * atoOut.second->capacity/100.0<<" amount interest paid to: "<<fid2<<" out of "<<nodeNum<<endl;
	// 		ir_pay[fid2] += rate * atoOut.second->capacity;
	// 		// ir_pay[fid2] += (int) atoOut.second->capacity/(10-atoOut.second->interest_rate);
	// 		// cout<<ir_pay[fid2]<<" will be paid to "<<fid2<<endl;
	// 	}
	// }
	// // for (int i=0; i< nodeSize; i++){
	// // 	checkCollateral(i);
	// // }	
	// // cout<< "added up "<<fid1<<endl;
	// for (int i=0; i< nodeSize; i++){
	// 	if (this->nodes[i]->defaulted || this->nodes[i]->isMarket){
	// 		// cout<<"node "<<i<< " is market or is already defaulted"<<endl;
	// 		continue;
	// 	}		
	// 	// cout<<ir_pay[i]<<" will be paid to "<<i<<endl;		
	// 	if(ir_pay[i]>0){
	// 		// cout<<"ir to pay by "<< fid1 << " to "<< i<< " is "<< ir_pay[i]<<" rounded to " << round(ir_pay[i]*precision)/precision<<endl;
	// 		// this->updateNodeDegrees();
	// 		// cout<<"updated degree"<<endl;
	// 		string mode = this->nodes[fid1]->routePreference;

	// 		// this->print();

	// 		// fid1 = 2; fid2 = 0;
	// 		// cout << "from " << fid1 << ", to " << fid2 << " " << request << " " << mode;

	// 		CplexConverter converter;
	// 		converter.constructCplex(this, this->nodes[fid1], this->nodes[i], round(ir_pay[i]*precision)/precision, transactionCounter);
	// 		transactionCounter += 1;
	// 		// converter.printInput();
			
	// 		LpSolver lpSolver;

			
	// 		if (lpSolver.solveLpProblem(converter, mode)){
	// 			// converter.printResult();
	// 			double check = converter.copyBack();
	// 			if (round(precision*check)/precision != round(precision*converter.request)/precision){
	// 				cout<< "interest amt check failed"<<endl;
	// 				cout<< "owed "<<ir_pay[i]<< " paid "<<check<<endl;
	// 				// cout<<"default"<<endl;
	// 				// this->nodes[fid1]->makeDefault();
	// 				// unwind(fid1);					
	// 				// return 1;
	// 			}				
	// 			// cout << " success " << endl;
	// 		}
	// 		// cout << " fail " << endl;
	// 		else{
	// 			// cout<<"default"<<endl;
	// 			// this->nodes[fid1]->print();	

	// 			// this->nodes[i]->print();	
	// 			// this->print();
	// 			this->nodes[fid1]->makeDefault();
	// 			unwind(fid1);				
	// 			// cout<<"deleted"<<endl;
	// 			// this->nodes[fid1]->print();
	// 			// cout<<"printed"<<endl;
	// 			return 1;
	// 		}
	// 	}
	// }
	return 0;
}
