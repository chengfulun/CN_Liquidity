#include "CN_CplexConverter.h"
#include <climits>
#include <map>
#include <algorithm>

extern CredNetConstants credNetConstants;

// 1 for sending, 0 for receiving
double CplexConverter::valueCluster(int val, int type){
	// remember that clusters must include 0 and max valuation
	if (type == 1){
	// sending needs to assume a valuation cluster lower than actual value
	// goal is to prevent all transactions that are not monotone
		// even if some monotone ones are also prevented
		double last_cluster = credNetConstants.totalValues[0];
		double cluster = 0;
		for (int i = 1; i < credNetConstants.totalValues.size(); ++i){
			cluster = credNetConstants.totalValues[i];
			if (val < cluster){
				return last_cluster;
			}
			last_cluster = cluster;
		}
		return cluster;
	}
	// receiving needs to assume a valuation cluster greater than value
	double cluster = INT_MIN;
	for (int i; i < credNetConstants.totalValues.size(); ++i){
		double cluster = credNetConstants.totalValues[i];
		if (val <= cluster){
			return cluster;
		}
	}	
}

void CplexConverter::constructCplex(Graph* g, Node* s, Node* t, double req, int transSeqNumT, bool verb){
	graph = g;
	src = s;
	dest = t;
	request = req;
	transSeqNum = transSeqNumT;
	this->verb = verb;

	int globalVarId = 0;
	for (int i = 0; i < g->atomicEdges.size(); i++){
		AtomicEdge* at = g->atomicEdges[i];

		if (at->capacity <= 0.0){
			continue;
		}
		// cout<<"atomic edge done"<<endl;
		bool leveraged = at->isLeveraged;
		bool defaultedTo = at->nodeTo->defaulted;
		bool defaultedFrom = at->nodeFrom->defaulted;
		int credId = at->singleCreditIndex;
		bool active = at->originEdge->singleCreditEdges[credId]->active;

		// bool intermediate = true;
		if (at->isDebt){
			// if(at->originEdge->nodeTo->nodeId==dest->nodeId ||at->originEdge->nodeTo->nodeId==(g->nodes.size()-1)){
			// 	intermediate = false;
			// }
			
			int tempId = globalVarId++;

			// for (int i; i < credNetConstants.totalValues.size(); ++i){
			// 	int ix = credNetConstants.totalValues.size() - 1 - i;
			// 	double val = credNetConstants.totalValues[ix];
			// 	double receiving_ir = val - at->collateral_value;
			// 	double sending_ir = (val +  (1-at->CR()) * at->default_rate) / (1-at->default_rate);
			// 	if (sending_ir >= receiving_ir and sending_ir ){
			// 		Varoable v(tempId,at->atomicEdgeId,val,)
			// 	}
			
			double receiving_val = this->valueCluster(at->getReceiveValue(at->interest_rate),0);
			double sending_val = this->valueCluster(at->getSendValue(at->interest_rate),1);

			// if (receiving_val <= sending_val){
			// 	Variable v(tempId, at->atomicEdgeId, receiving_val, sending_val, at->interest_rate);
			// 	variables.push_back(v);				
			// }
			Variable v(tempId, at->atomicEdgeId, receiving_val, sending_val, at->interest_rate);
			variables.push_back(v);

			if (defaultedTo){
				// || (intermediate &&(at->interest_rate==credNetConstants.totalIrs[1]) )){
				capacities.push_back(0.0);
			}
			else{
				capacities.push_back(max(0.0,at->capacity));
			}
			atomicIdToVarIdDict[at->atomicEdgeId].push_back(tempId);
		// cout<<"atov done"<<endl;

			varIdToAtomicIdDict[tempId] = at->atomicEdgeId;
					// cout<<"vtoa done"<<endl;


		}

		else {
			// if(at->originEdge->nodeFrom->nodeId==dest->nodeId || at->originEdge->nodeTo->nodeId==(g->nodes.size()-1)){
			// 	intermediate = false;
			// }
			double receiving_val = this->valueCluster(at->getReceiveValue(at->interest_rate),0);
			double sending_val = this->valueCluster(at->getSendValue(at->interest_rate),1);

			double diff = receiving_val - sending_val;

			for (int i = 0; i < credNetConstants.totalValues.size(); ++i){
				// routing value must exceed receiver's cluster
				if (credNetConstants.totalValues[i] > receiving_val){

					int tempId = globalVarId++;
					// cout << "adding credit: " << tempId << " " << at->atomicEdgeId << " " << at->interest_rate << endl;
					double paid_ir = (credNetConstants.totalValues[i] + (1-at->CR()) * at->default_rate) / (1-at->default_rate);
					// at->interest_rate + max(0,credNetConstants.totalValues[i] - at->getReceiveValue());
					double sent_value = this->valueCluster(paid_ir + at->collateral_value,1);
					// this->valueCluster(at->getSendValue() + at->interest_rate - paid_ir,1);
					Variable v(tempId, at->atomicEdgeId, credNetConstants.totalValues[i],sent_value , paid_ir);
					variables.push_back(v);
					if (defaultedFrom || leveraged || not active){
						// || (intermediate &&(at->interest_rate==credNetConstants.totalIrs[1]) )){
						capacities.push_back(0.0);
					}
					else{
						capacities.push_back(max(0.0,at->capacity));
					}
					atomicIdToVarIdDict[at->atomicEdgeId].push_back(tempId);
							// cout<<"atov done"<<endl;

					varIdToAtomicIdDict[tempId] = at->atomicEdgeId;
										// cout<<"vtoa done"<<endl;


				}
			}
		}
	}
}

void CplexConverter::printInput(){
	cout << "printing input variables... " << this->transSeqNum << endl;
	for (int i = 0; i < variables.size(); ++i){
		cout << "Var ID " << variables[i].varId 
		<< ", Atomic Edge ID " << variables[i].atomicEdgeId
		<< ", IR " << variables[i].interest_rate 
		<< ", range 0 - " << capacities[i] << endl;
	}
}

void CplexConverter::printResult(){
	for (int i = 0; i < results.size(); ++i){
		cout << results[i] << " ";
	}
	cout << endl;
}

double CplexConverter::getRouted(){
	double routed = 0;
	int srcId = this->src->nodeId;
	for (int i = 0; i < variables.size(); ++i){
		if (results[i] == 0){
			continue;
		}
		int id = variables[i].atomicEdgeId;		
		bool isDebt = this->graph->atomicEdges[id]->isDebt;

		int Id;
		int toId;
		if (isDebt){
			Id = this->graph->atomicEdges[id]->nodeFrom->nodeId;
			toId = this->graph->atomicEdges[id]->nodeTo->nodeId;
		}
		else{
			Id = this->graph->atomicEdges[id]->nodeTo->nodeId;
			toId = this->graph->atomicEdges[id]->nodeFrom->nodeId;			
		}
		if(Id == srcId){
			routed += results[i];
		}
	}
	return routed;
}

double CplexConverter::getAssetCost(){
	double capcost = 0;
	int srcId = this->src->nodeId;
	for (int i = 0; i < variables.size(); ++i){
		if (results[i] == 0){
			continue;
		}
		int id = variables[i].atomicEdgeId;	
		int ir = variables[i].interest_rate;			
		bool isDebt = this->graph->atomicEdges[id]->isDebt;

		int Id;
		int toId;
		if (isDebt){
			Id = this->graph->atomicEdges[id]->nodeFrom->nodeId;
			toId = this->graph->atomicEdges[id]->nodeTo->nodeId;
		}
		else{
			Id = this->graph->atomicEdges[id]->nodeTo->nodeId;
			toId = this->graph->atomicEdges[id]->nodeFrom->nodeId;			
		}
		if(Id == srcId){
			capcost += results[i] * this->graph->atomicEdges[id]->getSendValue(ir);
			// variables[i].sending_val;
		}
	}
	return capcost;
}

double CplexConverter::getDebtCost(){
	map<int, vector<double>> credits_used;

	double capcost = 0;
	int srcId = this->src->nodeId;
	int destId = this->dest->nodeId;
	for (int i = 0; i < variables.size(); ++i){
		int id = variables[i].atomicEdgeId;
		if (results[i] == 0){
			continue;
		}
		bool isDebt = this->graph->atomicEdges[id]->isDebt;

		int Id;
		int toId;
		if (isDebt){
			Id = this->graph->atomicEdges[id]->nodeFrom->nodeId;
			toId = this->graph->atomicEdges[id]->nodeTo->nodeId;
		}
		else{
			Id = this->graph->atomicEdges[id]->nodeTo->nodeId;
			toId = this->graph->atomicEdges[id]->nodeFrom->nodeId;			
		}
		if(toId != destId && not isDebt){
			if(credits_used.find(toId)==credits_used.end()){
				vector<double> new_entry;
				new_entry.push_back(results[i]);
				credits_used[toId] = new_entry;
			}
			else{
				credits_used[toId].push_back(results[i]);				
			}
		}
	}

	for (int i = 0; i < variables.size(); ++i){
		int Id;
		int toId;
		int id = variables[i].atomicEdgeId;
		bool isDebt = this->graph->atomicEdges[id]->isDebt;
		int ir = variables[i].interest_rate;			

		if (results[i] == 0){
			continue;
		}


		if (isDebt){
			Id = this->graph->atomicEdges[id]->nodeFrom->nodeId;
			toId = this->graph->atomicEdges[id]->nodeTo->nodeId;
		}
		else{
			Id = this->graph->atomicEdges[id]->nodeTo->nodeId;
			toId = this->graph->atomicEdges[id]->nodeFrom->nodeId;			
		}

		if(credits_used.find(Id)==credits_used.end()){
			continue;
		}

		if(Id != srcId){

			// std::vector<double> this_credits = credits_used[Id];

			std::vector<double>::iterator position = std::find(credits_used[Id].begin(), credits_used[Id].end(), results[i]);
			if (position != credits_used[Id].end()){
				capcost += results[i] * this->graph->atomicEdges[id]->getSendValue(ir);			
		    	credits_used[Id].erase(position);			
			} // == myVector.end() means the element was not found

			// for (auto it = this_credits.begin();it != this_credits.end(); it++){
			// 	if(results[i] == *it){
			// 		// std::cout<<it<<endl;
			// 		capcost += results[i] * this->graph->atomicEdges[id]->getSendValue(ir);
			// 		// capcost += results[i] * variables[i].value_paying;
			// 		it.erase(it--);
			// 		break;
			// 	}
			// }

		}
	}
	return capcost;
} 

double CplexConverter::copyBack(){
	double routed = 0;
	int srcId = this->src->nodeId;
	int destId = this->dest->nodeId;	
	// cout<<"copying source "<<srcId<<" dest "<<destId<<endl;
	for (int i = 0; i < variables.size(); ++i){
		int id = variables[i].atomicEdgeId;
		if (results[i] == 0){
			continue;
		}
		bool isDebt = this->graph->atomicEdges[id]->isDebt;

		int Id;
		int toId;
		if (isDebt){
			Id = this->graph->atomicEdges[id]->nodeFrom->nodeId;
			toId = this->graph->atomicEdges[id]->nodeTo->nodeId;
		}
		else{
			Id = this->graph->atomicEdges[id]->nodeTo->nodeId;
			toId = this->graph->atomicEdges[id]->nodeFrom->nodeId;			
		}
		if(Id == srcId){
			routed += results[i];
		}
		if(verb){
			if(Id != this->graph->marketId && toId!=this->graph->marketId){
				cout<<"routed "<< results[i]<<" credit to "<<toId<< " from "<<Id<< " at IR "<<variables[i].interest_rate<<endl;
			}
			if(toId == this->graph->marketId){
				cout<<"routed "<< results[i]<<" collateral from "<<Id<<endl;
			}
			if(Id == this->graph->marketId){
				cout<<"routed "<< results[i]<<" collateral to "<<toId<<endl;
			}	
		}	
		// cout<<"copied ir: "<<variables[i].print()<<endl;
		this->graph->atomicEdges[id]
			->route(results[i], variables[i].interest_rate, this->graph, this->transSeqNum);
	}
	// cout<<"copying done"<<endl;
	return routed;
}