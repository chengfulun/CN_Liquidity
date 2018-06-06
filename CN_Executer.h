#ifndef CN_Executer
#define CN_Executer

#include "CN_Node.h"
#include "CN_Graph.h"
#include <list>

class Executer{
public:
	static void execute(vector<Edge*>& path, Node* src, Node* dest);
};


#endif