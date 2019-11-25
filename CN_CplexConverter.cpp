#include "CN_CplexConverter.h"

extern CredNetConstants credNetConstants;

// 1 for sending, 0 for receiving
double CplexConverter::valueCluster(int val, int type){
	if (type == 0){
	// receiving needs to assume a valuation cluster lower than actual value
		double last_cluster = INT_MIN;
		double cluster = 0;
		for (int i; i < credNetConstants.totalValues.size(); ++i){
			double cluster = credNetConstants.totalValues[i];
			if (val <= cluster){
				return last_cluster;
			}
			last_cluster = cluster;
		}
	}
	// sending needs to assume a valuation cluster greater than value
	double cluster = INT_MIN;
	for (int i; i < credNetConst`ants.totalValues.size(); ++i){
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
		if (at->capacity <= 0.0){
			continue;
		}		
		AtomicEdge* at = g->atomicEdges[i];
		// cout<<"atomic edge done"<<endl;
		bool leveraged = at->isLeveraged;
		bool defaultedTo = at->nodeTo->defaulted;
		bool defaultedFrom = at->nodeFrom->defaulted;
		bool active = at->singleCreditEdges[credId]->active;

		// bool intermediate = true;
		if (at->isDebt){
			// if(at->originEdge->nodeTo->nodeId==dest->nodeId ||at->originEdge->nodeTo->nodeId==(g->nodes.size()-1)){
			// 	intermediate = false;
			// }
			int tempId = globalVarId++;
			double receiving_val = this->valueCluster(at->getReceiveValue(),0);
			double sending_val = this->valueCluster(at->getSendValue(),1);

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
			double receiving_val = this->valueCluster(at->getReceiveValue(),0);
			double sending_val = this->valueCluster(at->getSendValue(),1);

			double diff = at->getReceiveValue() - at->getSendValue();

			for (int i = 0; i < credNetConstants.totalValues.size(); ++i){
				if (credNetConstants.totalValues[i] >= receiving_val){

					int tempId = globalVarId++;
					// cout << "adding credit: " << tempId << " " << at->atomicEdgeId << " " << at->interest_rate << endl;
					double paid_ir = at->interest_rate + max(0,credNetConstants.totalValues[i] - at->getReceiveValue());
					double sent_value = this->valueCluster(at->getSendValue() + at->interest_rate - paid_ir,1);
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
			capcost += results[i] * variables[i].sending_val;
		}
	}
	return capcost;
}

// modify this for when payments are done to return debt or deposits
double CplexConverter::getDebtCost(){
	double capcost = 0;
	int srcId = this->src->nodeId;
	for (int i = 0; i < variables.size(); ++i){
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
		if(Id != srcId){
			capcost += results[i] * variables[i].sending_val;
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