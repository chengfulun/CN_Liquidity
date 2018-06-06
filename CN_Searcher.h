#ifndef CN_Searcher
#define CN_Searcher

#include "CN_Node.h"
#include "CN_Graph.h"

#include <list>
#include <vector>

class Searcher{
public:

	static bool bfsIRConstraint(double ir, Graph* graph, 
		Node* src, Node* dest, vector<AtomicEdge*>& path);
	
};

#endif