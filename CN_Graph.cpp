
#include "CN_Graph.h"
#include "CN_Constants.h"
#include <iostream>
#include <random>
#include <algorithm>
#include <cfloat>
#include <limits>
#include <string>
#include <sstream>
#include <iostream>

using namespace std;

extern CredNetConstants credNetConstants;

/////////////////////////////////////////////////////////////////////////
/* Graph basics */
/////////////////////////////////////////////////////////////////////////
Graph::Graph(){
	this->atomicGlobalId = 0;
}

Graph::Graph(int nodeNumT, int marketId){
	this->nodeNum = nodeNumT;
	this->atomicGlobalId = 0;
	for (int i = 0; i < nodeNum; ++i){
		Node* temp = new Node(i);
		std::pair<int, Node*> tempPair;
		tempPair.first = i;
		tempPair.second = temp;
		nodes.insert(tempPair);
		if (i == marketId){
			nodes[i]->makeMarket();
		}
	}
	this->marketId = marketId;
}

// static void deepCopyHelper(Graph* newGraph, Graph& oldGraph){

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

// Graph::Graph(Graph &graphT){
// 	deepCopyHelper(this, graphT);
// }

// Graph& Graph::operator=(Graph &graphT){
// 	deepCopyHelper(this, graphT);
// 	return *this;
// }

// Graph::~Graph(){
// 	for (& it : nodes){
// 		delete it.second;
// 	}
// }

// void Graph::deleteEdges(){
// 	for (auto& it : edges){
// 		if (!it->nodeFrom->isMarket){
// 			delete it;			
// 		}
// 	}
// 	for (auto& it : atomicEdges){
// 		if (!it->nodeFrom->isMarket){
// 			delete it.second;
// 		}
// 	}
// }

void Graph::updateNodeDegrees(){
	for (auto& it : this->nodes){
		it.second->updateDegree();
	}
}

void Graph::modifyCredit(int source,int dest,double credit){
	// double sumC = 0.0;
	// cout<<" amt to set credit "<<credit<<endl;
	for (auto& e : this->nodes[source]->edge_out){
		bool done = false;
		if(done){
			break;
		}
		if(e.second->nodeTo->nodeId == dest){
			for (auto& ce : e.second->singleCreditEdges){
				if(not done){
					ce->setCredit(credit);
					done = true;
				}
				else{
					ce->setCredit(0.0);
				}

			// 	double sumNext = sumC + ce->credit_max;
			// 	if(sumNext > credit){
			// 		ce->setCredit(max(0.0,credit - sumC));
			// 	}
			// 	sumC += ce->credit_max;
			// }
			// if(sumC < credit){
			// 	e.second->singleCreditEdges[0]->setCredit(e.second->singleCreditEdges[0]->credit_max + credit - sumC);
			// }
			}
		}
	}
}

void Graph::addCredit(int source,int dest,double credit){
	// double sumC = 0.0;
	// cout<<" amt to set credit "<<credit<<endl;
	for (auto& e : this->nodes[source]->edge_out){
		bool done = false;
		if(done){
			break;
		}
		if(e.second->nodeTo->nodeId == dest){
			for (auto& ce : e.second->singleCreditEdges){
				if(not done){
					ce->setCredit(credit + ce->credit_max);
					done = true;
				}
				else{
					ce->setCredit(0.0);
				}

			// 	double sumNext = sumC + ce->credit_max;
			// 	if(sumNext > credit){
			// 		ce->setCredit(max(0.0,credit - sumC));
			// 	}
			// 	sumC += ce->credit_max;
			// }
			// if(sumC < credit){
			// 	e.second->singleCreditEdges[0]->setCredit(e.second->singleCreditEdges[0]->credit_max + credit - sumC);
			// }
			}
		}
	}
}

void Graph::addMultiEdge(Node* nodeFrom, Node* nodeTo, 
	double credit_ir, double debt_ir, double currDebt, double cap, double cr){

	// use diff ir to create ious
	Edge* e;
	if (nodeTo->edge_in.find(nodeFrom->getNodeId()) == nodeTo->edge_in.end()){
		e = new Edge(nodeFrom, nodeTo);
		nodeTo->edge_in[nodeFrom->getNodeId()] = e;
		nodeFrom->edge_out[nodeTo->getNodeId()] = e;
		this->edges.push_back(e);
	} else {
		e = nodeTo->edge_in[nodeFrom->getNodeId()];
	}

	SingleCreditEdge* temp = e->addSingleCreditEdge(credit_ir, cap, this->atomicGlobalId, this->atomicEdges, cr);
	if (currDebt != 0){
		e->routeAtomicEdge(temp->credit_remain, currDebt, debt_ir, atomicGlobalId, this->atomicEdges, 0);
	}

}

// void Graph::updateCreditEdgeCapacities(ix, vector<double> caps, ){

// }

void Graph::print(){
	for (auto& it : nodes){
		it.second->print();
	}
}

void Graph::printAtomicEdges(){
	cout << "--atomic edges of graph: ---------------------------------------" << endl;
	for (auto &it : this->atomicEdges){
		it.second->print();
	}
	cout << endl;

}

void Graph::printAvgAtomicIouEdges(){
	int cnt = 0;
	double totalIrCap = 0;
	double totalCap = 0;

	for (auto &it : this->atomicEdges){
		if (!it.second->isDebt || it.second->capacity == 0){
			continue;
		}
		totalIrCap += it.second->capacity * it.second->interest_rate;
		totalCap += it.second->capacity;
		cnt ++;
	}
	cout << "   " << totalIrCap / totalCap << endl;
}

void Graph::printAtomicIouEdges(ofstream & fout){
	fout << "--All IOU Edges ----------------------------------------------" << endl;
	fout << "nodeFrom \t nodeTo \t amount \t ir" << endl;
	for (auto &it : this->atomicEdges){
		if (!it.second->isDebt || it.second->capacity == 0){
			continue;
		}
		fout << it.second->nodeFrom->nodeId << " \t "
			<< it.second->nodeTo->nodeId << " \t "
			<< it.second->capacity << " \t "
			<< it.second->interest_rate << endl;
	}
	fout << endl;
}




/////////////////////////////////////////////////////////////////////////
/* Generate Initial Network */
/////////////////////////////////////////////////////////////////////////
void Graph::genMarket0Graph(double deposit, double shock, double wealth, int marketId, bool setT){

	for (int i = 0; i < nodeNum; ++i){
		Node* temp = new Node(i);
		// Node* temp = new Node(i);

		std::pair<int, Node*> tempPair;
		tempPair.first = i;
		tempPair.second = temp;
		nodes.insert(tempPair);
		if (i == marketId){
			nodes.find(i)->second->makeMarket();
		}		
	}

	this->marketId = marketId;	
	this->expected_deposit = deposit;
	this->deposit_shock = shock;
	this->FFR = FFR;
	this->init_wealth = wealth;
		// default_random_engine generator;
	// uniform_real_distribution<double> distribution(0.0, 1.0);
	// cout<<"nodes done"<<endl;
	for (int i = 0; i < nodeNum; i++){
		if(setT){
			double theta_draw = (credNetConstants.uniformDoubleDistribution(
								credNetConstants.globalGenerator));
			// double deposit_shock = (credNetConstants.uniformDoubleDistribution(
			// 					credNetConstants.globalGenerator))*deposit/shock;

			this->nodes[i]->theta = theta_draw;
		}// cout<<this->nodes[i]->theta<<endl;
		this->nodes[i]->deposits = this->expected_deposit;
		// if (rand()%2 == 1) {			
		// 	this->nodes[i]->deposits = this->nodes[i]->deposits + 2 * deposit_shock;
		// }

		// this->nodes[i]->assets[0] = 
		// this->nodes[i]->invest();
	}
	for (int i = 0; i < nodeNum; i++){
		for(int j = 0; j < nodeNum; j++){


			if (i != marketId){
				if (j == marketId){
					double currency = wealth + this->nodes[i]->deposits;

					this->addMultiEdge(nodes.find(i)->second, nodes.find(marketId)->second,0.0,0.0,currency,DBL_MAX,0.0);
				}
				else if(j != i){
					this->addMultiEdge(nodes.find(i)->second, nodes.find(j)->second,FFR,0.0,0.0,0.0,0.0);
				}
				// this->returns.push_back(FFR);
				// this->volatilities.push_back(0.2);
				// this->wealths.push_back(wealth);
			}
		
		// for (int j = 0; j < nodeNum; j++){
		// 	if (i != marketId && j != marketId && i != j){

		// 		// double ir = (credNetConstants.uniformIntDistribution(
		// 		// 	credNetConstants.globalGenerator) % numIR);
		// 		// cout<<"random: "<<ir<<endl;	
		// 			// this->addUnitEdge(nodes.find(i)->second, nodes.find(j)->second, ir, rand()%2);
		// 		double cap = this->nodes[i]->getCash()/(nodeNum - 1);
		// 		this->addMultiEdge(nodes.find(i)->second, nodes.find(j)->second, FFR, 0.0, 0, cap, CR);
		// 		// this->addUnitEdge(nodes.find(j)->second, nodes.find(i)->second, ir, rand()%2);
		// 	}

		// }
		}
	}
	// cout<<"edges done"<<endl;


}

void Graph::genTest0Graph(double threshold, int numIR, double cap, double maxCR,double wealth, int marketId){

	for (int i = 0; i < nodeNum; ++i){
		Node* temp = new Node(i);
		std::pair<int, Node*> tempPair;
		tempPair.first = i;
		tempPair.second = temp;
		nodes.insert(tempPair);
		if (i == marketId){
			nodes.find(i)->second->makeMarket();
		}		
	}

	this->marketId = marketId;	

	// default_random_engine generator;
	// uniform_real_distribution<double> distribution(0.0, 1.0);
	// cout<<"nodes done"<<endl;
	for (int i = 0; i < nodeNum; i++){
		for (int j = i+1; j < nodeNum; j++){
			if (i != marketId && j != marketId){
				double num = credNetConstants.uniformDoubleDistribution(
					credNetConstants.globalGenerator);

				double ir = (credNetConstants.uniformIntDistribution(
					credNetConstants.globalGenerator) % numIR);
				// cout<<"random: "<<ir<<endl;
				double cr = (credNetConstants.uniformDoubleDistribution(
					credNetConstants.globalGenerator) * maxCR);
				
				if (num < threshold){
					if (rand()%2 == 1) {
						// this->addUnitEdge(nodes.find(i)->second, nodes.find(j)->second, ir, rand()%2);
						this->addMultiEdge(nodes.find(i)->second, nodes.find(j)->second, ir, 0.0, 0, cap, cr);
					} else {
						// this->addUnitEdge(nodes.find(j)->second, nodes.find(i)->second, ir, rand()%2);
						this->addMultiEdge(nodes.find(j)->second, nodes.find(i)->second, ir, 0.0, 0, cap, cr);
					}
				}				
			}

		}
	}


	for (int j = 0; j<nodeNum; j++){
		if (j != marketId){
			this->addMultiEdge(nodes.find(j)->second, nodes.find(marketId)->second,0.0,0.0,wealth,DBL_MAX,0.0);
		}
		else{
			this->nodes[j]->makeMarket();
		}
	}

	// cout<<"edges done"<<endl;


}


// void Graph::genTest1Graph(double threshold, int numIR, int cap){

// 	for (int i = 0; i < nodeNum; ++i){
// 		Node* temp = new Node(i);
// 		std::pair<int, Node*> tempPair;
// 		tempPair.first = i;
// 		tempPair.second = temp;
// 		nodes.insert(tempPair);
// 	}

// 	// default_random_engine generator;
// 	// uniform_real_distribution<double> distribution(0.0, 1.0);
// 	for (int i = 0; i < nodeNum; i++){

// 		int count = threshold / 2 * nodeNum + 1;
// 		while(count > 0){
			
// 			double ir = (credNetConstants.uniformIntDistribution(
// 				credNetConstants.gloabalGenerator) % numIR + 1);

// 			double id = (credNetConstants.uniformIntDistribution(
// 				credNetConstants.gloabalGenerator) % nodeNum);
// 			if (id == i){
// 				continue;
// 			}
// 			this->addMultiEdge(nodes.find(id)->second, nodes.find(i)->second, ir, 0.0, 0, cap);
// 			count--;
// 		}

// 	}
// }


// void Graph::generateTestGraph2(){

// 	nodeNum = 8;
// 	for (int i = 0; i < 8; ++i){
// 		Node* temp = new Node(i);
// 		std::pair<int, Node*> tempPair;
// 		tempPair.first = i;
// 		tempPair.second = temp;
// 		nodes.insert(tempPair);
// 	}
// 	this->addMultiEdge(nodes.find(0)->second, nodes.find(1)->second, 0.2, 0.4, 2, 2);

// 	this->addMultiEdge(nodes.find(0)->second, nodes.find(1)->second, 0.3, 0.3, 3, 3);

// 	this->addMultiEdge(nodes.find(2)->second, nodes.find(1)->second, 0.3, 0.0, 0, 2);

// 	this->addMultiEdge(nodes.find(2)->second, nodes.find(4)->second, 0.2, 0.2, 4, 4);

// 	this->addMultiEdge(nodes.find(3)->second, nodes.find(0)->second, 0.2, 0.0, 0, 2);

// 	this->addMultiEdge(nodes.find(4)->second, nodes.find(3)->second, 0.4, 0.0, 0, 1);

// 	this->addMultiEdge(nodes.find(5)->second, nodes.find(4)->second, 0.2, 0.0, 0, 3);

// 	this->addMultiEdge(nodes.find(5)->second, nodes.find(2)->second, 0.4, 0.0, 0, 2);

// 	this->addMultiEdge(nodes.find(6)->second, nodes.find(1)->second, 0.2, 0.0, 0, 2);

// 	this->addMultiEdge(nodes.find(7)->second, nodes.find(6)->second, 0.3, 0.0, 0, 1);

// 	this->addMultiEdge(nodes.find(4)->second, nodes.find(7)->second, 0.2, 0.0, 0, 2);

// }


// void Graph::generateTestGraph3(){

// 	nodeNum = 4;
// 	for (int i = 0; i < nodeNum; ++i){
// 		nodes[i] = new Node(i);
// 	}

// 	this->addMultiEdge(nodes.find(0)->second, nodes.find(1)->second, 0.02, 0.0, 0, 2);

// 	this->addMultiEdge(nodes.find(0)->second, nodes.find(2)->second, 0.03, 0.0, 0, 2);

// 	this->addMultiEdge(nodes.find(1)->second, nodes.find(0)->second, 0.02, 0.0, 0, 2);

// 	this->addMultiEdge(nodes.find(1)->second, nodes.find(2)->second, 0.01, 0.0, 0, 2);

// 	this->addMultiEdge(nodes.find(1)->second, nodes.find(3)->second, 0.01, 0.0, 0, 2);

// 	this->addMultiEdge(nodes.find(2)->second, nodes.find(0)->second, 0.03, 0.0, 0, 2);

// 	this->addMultiEdge(nodes.find(3)->second, nodes.find(2)->second, 0.02, 0.0, 0, 2);

// }



void Graph::setRoutePreference(vector<string> &v){
     
	for (int i = 0; i < nodeNum; i++) {
		nodes[i]->routePreference = v[i];
	}

	return;
}


void Graph::setThetas(vector<string> &v){
	for (int i = 0; i < nodeNum - 1; i++) {
		stringstream ss(v[i]);
		vector<string> res;

		while(ss.good() ){
			string substr;
			getline(ss,substr, '_');
			res.push_back(substr);
		}
		nodes[i]->theta = stod(res[0]);
		nodes[i]->CR_base = stod(res[1]);
	}

	return;
}

void Graph::printThetas(){
	for (int i = 0; i < nodeNum - 1; i++){
		std::cout<<"id: "<<i<<" theta: "<<nodes[i]->theta<<" CR: "<<nodes[i]->CR_base<<endl;
	}
}
// void Graph::setThetas(vector<double> &v){
     
// 	for (int i = 0; i < nodeNum - 1; i++) {
// 		nodes[i]->theta = v[i];
// 	}

// 	return;
// }

//////////////////////////////////////////////////////
// route on atomic edge
//////////////////////////////////////////////////////
void helpRouteOnAtomicEdge(double current, double interest_rate, AtomicEdge* a, Graph* g, int transSeqNum){
	a->originEdge->routeAtomicEdge(a, current, interest_rate, g->atomicGlobalId, g->atomicEdges, transSeqNum);
}
