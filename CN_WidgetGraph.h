#ifndef CN_WidgetGraph
#define CN_WidgetGraph

#include "CN_Node.h"
#include "CN_Graph.h"
#include <unordered_map>
#include <iostream>
#include <string>
#include <vector>
#include <set>


using namespace std;

enum WidgetNodeType{
	CREDIT_IN_NODE, // Out Widget Edge
	DEBT_OUT_NODE, // Out Widget Edge
	CREDIT_OUT_NODE, // In Widget Edge
	DEBT_IN_NODE, // In Widget Edge
	SUPER_SOURCE,
	SUPER_DEST,
};

enum WidgetEdgeType{
	CREDIT_WIDGET_EDGE,
	DEBT_WIDGET_EDGE,
	INNER_WIDGET_EDGE,
};

static std::string helper(WidgetNodeType type){
	switch(type){
		case CREDIT_OUT_NODE:
			return "CREDIT_OUT_NODE";
		case CREDIT_IN_NODE:
			return "CREDIT_IN_NODE";
		case DEBT_OUT_NODE:
			return "DEBT_OUT_NODE";
		case DEBT_IN_NODE:
			return "DEBT_IN_NODE";
	}
	cerr << "unknown widget node type";
	return "";
}

static std::string helper(WidgetEdgeType type){
	switch(type){
		case CREDIT_WIDGET_EDGE:
			return "CREDIT_WIDGET_EDGE";
		case DEBT_WIDGET_EDGE:
			return "DEBT_WIDGET_EDGE";
		case INNER_WIDGET_EDGE:
			return "INNER_WIDGET_EDGE";
	}
	cerr << "unknown widget edge type";
	return "";
}

class WidgetNode;

struct WidgetEdge{

	AtomicEdge* originAtomicEdge;

	WidgetNode* nodeFrom;
	WidgetNode* nodeTo;
	double cap;
	double curr;
	double interest_rate;
	double interest_diff;

	WidgetEdgeType type;

	WidgetEdge(double ir, double ir_diff, int capT, WidgetNode* nodeFromT, 
		WidgetNode* nodeToT, WidgetEdgeType typeT, AtomicEdge* a)
		: curr(0), cap(capT), interest_rate(ir), nodeTo(nodeToT)
		, nodeFrom(nodeFromT), interest_diff(ir_diff), type(typeT), originAtomicEdge(a) {}
};

class WidgetNode{
public:
	int globalNodeId;
	WidgetNodeType type;
	Node* originNode;
	int targetNodeId;
	AtomicEdge* originAtomicEdge;

	unordered_map<int, WidgetEdge*> edge_out;
	unordered_map<int, WidgetEdge*> edge_in;

	WidgetNode(WidgetNodeType t, Node* o, int id, int targetNodeIdT, AtomicEdge* a)
		: type(t), originNode(o)
		, globalNodeId(id), targetNodeId(targetNodeIdT) 
		, originAtomicEdge(a) {}

	// ~WidgetNode(){
	// 	for (auto it : edge_out){
	// 		delete it.second;
	// 	}
	// }

	int getGlobalId(){
		return this->globalNodeId;
	}

	void print(){
		cout << "Global ID: " << this->globalNodeId << endl;
		cout << "origin node: " << originNode->getNodeId() << " " 
			<< helper(type) << ", target node: " << targetNodeId << endl;

		for (auto it : edge_in){
			cout << "From widget node " << it.second->nodeFrom->globalNodeId 
				<< " which belongs to " <<  it.second->nodeFrom->originNode->getNodeId() 
				<< " Type: " << helper(it.second->type) 
				<< endl
				<< "interest rate: " << it.second->interest_rate 
				<< " interest diff: " << it.second->interest_diff
				<< " capacity: " << it.second->cap 
				<< " current: " << it.second->curr << endl;
		}

		for (auto it : edge_out){
			cout << "To widget node " << it.second->nodeTo->globalNodeId 
				<< " which belongs to " <<  it.second->nodeTo->originNode->getNodeId() 
				<< " Type: " << helper(it.second->type) 
				<< endl
				<< "interest rate: " << it.second->interest_rate 
				<< " interest diff: " << it.second->interest_diff
				<< " capacity: " << it.second->cap 
				<< " current: " << it.second->curr << endl;
		}

		cout << endl;
	}
};

/**
 * all credit edges will be inversed, 
 * 
 */
class WidgetGraph{
private:
	Graph* originGraph;
public:
	vector<WidgetNode*> widgetNodes;
	set<WidgetEdge*> outerWidgetEdges;
	set<WidgetEdge*> widgetEdges;

	// super source and super sink
	Node* superSource;
	Node* superSink;
	double payment;

	void addEdge(WidgetNode* node1, WidgetNode* node2, 
		int capacity, double ir, double ir_diff, 
		WidgetEdgeType type, AtomicEdge* a, bool isOuter);
	// void constructWidget(Graph* graphT);
	void copyBack();

	WidgetGraph();
	// ~WidgetGraph();
	
	void print(){
		for (int i = 0; i < widgetNodes.size(); ++i){
			widgetNodes[i]->print();
		}
	}
	void setupSrcAndDest(Node*, Node*, double);
};

#endif