
#include "CN_Edge.h"
#include <iostream>

using namespace std;

Edge::Edge(Node* nodeFromT, Node* nodeToT){
	this->nodeFrom = nodeFromT;
	this->nodeTo = nodeToT;
}

Edge::Edge(const Edge& e, Node* nodeFromT, Node* nodeToT, 
	unordered_map<int, AtomicEdge*>& atomicMap){
	this->nodeFrom = nodeFromT;
	this->nodeTo = nodeToT;
	for (int i = 0; i < e.singleCreditEdges.size(); ++i){
		SingleCreditEdge* s = new SingleCreditEdge(*(e.singleCreditEdges[i]), 
			this, this->singleCreditEdges.size(), atomicMap, nodeFromT, nodeToT);
		this->singleCreditEdges.push_back(s);
	}
}

// Edge::~Edge(){
// 	for (int i = 0; i < singleCreditEdges.size(); ++i){
// 		delete singleCreditEdges[i];
// 	}
// }

SingleCreditEdge* Edge::addSingleCreditEdge(double interest_rate, double cap, 
	int& atomicGlobalId, unordered_map<int, AtomicEdge*>& atomicMap, double collateralRate){
	// cout<<"cred interest: "<<interest_rate<<endl;
		SingleCreditEdge* temp = new SingleCreditEdge(cap, 

	// SingleCreditEdge* temp = new SingleCreditEdge(cap, 
		interest_rate, atomicGlobalId, this, 
		this->singleCreditEdges.size(), atomicMap, nodeFrom, nodeTo, collateralRate);
	this->singleCreditEdges.push_back(temp);
	return temp;
	
}

void Edge::routeAtomicEdge(AtomicEdge* a, double flow, double interest_rate, 
	int& atomicGlobalId, unordered_map<int, AtomicEdge*>& atomicMap, int transSeqNum){
	// cout<<"debt ir: "<<interest_rate<<endl;
	if (a->isDebt){
		a->originEdge->singleCreditEdges[a->singleCreditIndex]
			->routeDebt(flow, a->interest_rate);
		a->nodeFrom->addModification(transSeqNum);
	} else {
		// cout<<"cred ir: "<<interest_rate<<endl;
		a->originEdge->singleCreditEdges[a->singleCreditIndex]
			->routeCredit(flow, interest_rate, atomicGlobalId, 
				this, a->singleCreditIndex, atomicMap, nodeFrom, nodeTo);
		a->nodeTo->addModification(transSeqNum);
	}

	this->edgeUsage.push_back(transSeqNum);

}

void Edge::print(){
	cout << "From Node " << nodeFrom->getNodeId() 
	<< " To Node " << nodeTo->getNodeId() << endl;

	for (int i = 0; i < singleCreditEdges.size(); ++i){
		singleCreditEdges[i]->print();
	}
	
}


void printEdge(Edge* e){
	e->print();
}

void printAtomicEdge(AtomicEdge* ae){
	ae->print();
}

double getNodeCurrBalance(Node* n){
	double temp = 0;
	for (auto& it : n->atomicEdge_in){
		if (it.second->isDebt){
			temp += it.second->capacity * it.second->interest_rate/100;
		}
	}
	// cout<<"debt in done"<<endl;
	for (auto& it : n->atomicEdge_out){
		if (it.second->isDebt){
			temp -= it.second->capacity * it.second->interest_rate/100;
		}
	}
	// cout<<"debt out done"<<endl;

	return temp;
}