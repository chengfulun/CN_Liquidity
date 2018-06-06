#include "CN_Searcher.h"

#define NODE_NUM_MAX 1000

#include <set>
#include <queue>
#include <vector>
#include <list>
#include <stack>
#include <unordered_map>
#include <iostream>

using namespace std;

static void deepCopyHelper(vector<Edge*>& a, vector<Edge*>& b){
	b.clear();
	for (int i = 0; i < a.size(); ++i){
		b.push_back(a[i]);
	}
}

struct BfsQueueItem{
	Node* front;
	double currIR;
	vector<Edge*> path;
	BfsQueueItem(Node* f, double ir, vector<Edge*>& pathT) : front(f), currIR(ir) {
		deepCopyHelper(pathT, path);
	}
};

bool Searcher::bfsIRConstraint(double ir, 
	Graph* graph, Node* src, Node* dest, vector<Edge*>& path){

	queue <BfsQueueItem*> tempQueue;
	unordered_map<int, double> tempVisited; // nodeId, max ir
	Node* front;


	// initialization
	for (auto edge : src->edge_in){

		if (edge.second->get_c_remain() == 0 ||
			ir < edge.second->get_interest_rate()){
			continue;
		}

		pair<int, double> visitedPair;
		visitedPair.first = edge.second->nodeFrom->getNodeId();
		visitedPair.second = ir;
		tempVisited.insert(visitedPair);

		vector<Edge*> tempPath;
		tempPath.push_back(edge.second);
		BfsQueueItem* queueItem = 
			new BfsQueueItem(edge.second->nodeFrom, ir, tempPath);
		tempQueue.push(queueItem);

	}
	for (auto edge : src->edge_out){

		if (edge.second->get_d_current() == 0 ||
			ir < edge.second->get_debt_interest_rate()){
			continue;
		}

		pair<int, double> visitedPair;
		visitedPair.first = edge.second->nodeTo->getNodeId();
		visitedPair.second = ir;
		tempVisited.insert(visitedPair);

		vector<Edge*> tempPath;
		tempPath.push_back(edge.second);
		BfsQueueItem* queueItem = 
			new BfsQueueItem(edge.second->nodeTo, edge.second->get_debt_interest_rate(), tempPath);
		tempQueue.push(queueItem);

	}
	
	// cout << "initialized " << tempQueue.size() << endl;

	while(tempQueue.size() != 0){

		BfsQueueItem* frontPair = tempQueue.front();
		double currIR = frontPair->currIR;
		front = frontPair->front;
		deepCopyHelper(frontPair->path, path);
		// cout << "path length: " << frontPair->path.size() << " " << path.size() << " " << endl;
		tempQueue.pop();
		delete frontPair;

		// cout << "front: " << front->getNodeId() << endl;
		// cout << "currIR: " << currIR << endl;
		// for (int i = 0; i < path.size(); ++i){
		// 	cout << path[i]->nodeFrom->getNodeId() << path[i]->nodeTo->getNodeId() << ", ";
		// }
		// cout << endl;

		if(dest == front){
			break;
		}

		// push to queue
		for(auto edge : front->edge_in){ 

			if (edge.second->get_c_remain() == 0){
				continue;
			}

			if(tempVisited.end() == tempVisited.find(edge.second->nodeFrom->getNodeId())){
				if (edge.second->get_interest_rate() > currIR){
					continue;
				}

				// cout << "front->edge " << edge.second->nodeTo->getNodeId() << " " << edge.second->nodeFrom->getNodeId() << endl;

				pair<int, double> visitedPair;
				visitedPair.first = edge.second->nodeFrom->getNodeId();
				visitedPair.second = currIR;
				tempVisited.insert(visitedPair);
			} else if (tempVisited.find(edge.second->nodeFrom->getNodeId())->second < currIR){

				// cout << "front->edge " << edge.second->nodeTo->getNodeId() << " " << edge.second->nodeFrom->getNodeId() << endl;

				tempVisited.find(edge.second->nodeFrom->getNodeId())->second = currIR;
			} else {
				continue;
			}

			vector<Edge*> tempPath;
			deepCopyHelper(path, tempPath);
			tempPath.push_back(edge.second);
			BfsQueueItem* queueItem = 
				new BfsQueueItem(edge.second->nodeFrom, currIR, tempPath);
			tempQueue.push(queueItem);

		}

		for(auto edge : front->edge_out){ 
			if (edge.second->get_d_current() == 0){
				continue;
			}

			if(tempVisited.end() == tempVisited.find(edge.second->nodeTo->getNodeId())){
				
				if (edge.second->get_debt_interest_rate() > currIR){
					continue;
				}

				// cout << "front->edge " << edge.second->nodeFrom->getNodeId() << " " << edge.second->nodeTo->getNodeId() << endl;

				pair<int, double> visitedPair;
				visitedPair.first = edge.second->nodeTo->getNodeId();
				visitedPair.second = edge.second->get_debt_interest_rate();
				tempVisited.insert(visitedPair);

			} else if (tempVisited.find(edge.second->nodeTo->getNodeId())->second < currIR){

				// cout << "front->edge " << edge.second->nodeFrom->getNodeId() << " " << edge.second->nodeTo->getNodeId() << endl;

				tempVisited.find(edge.second->nodeTo->getNodeId())->second = currIR;
			} else {
				continue;
			}

			vector<Edge*> tempPath;
			deepCopyHelper(path, tempPath);
			tempPath.push_back(edge.second);
			BfsQueueItem* queueItem = 
				new BfsQueueItem(edge.second->nodeTo, edge.second->get_debt_interest_rate(), tempPath);
			tempQueue.push(queueItem);

		}

// cout << "queue size: " << tempQueue.size() << endl;
	}

	while (!tempQueue.empty()) {
		delete tempQueue.front();
		tempQueue.pop();
	}

// cout << "path length: " << path.size() << endl;
	// cout << "front" << front->getNodeId() << endl;

	if(front != dest){
		path.clear(); 
		return false; 
	}

	return true;
}

