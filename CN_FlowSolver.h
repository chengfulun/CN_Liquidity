#ifndef CN_CplexSolver
#define CN_CplexSolver

#include "CN_WidgetGraph.h"
#include <vector>

class CplexConverter{

public:
	// widget graph
	int widgetNodeNum;
	int widgetEdgeNum;
	vector<WidgetEdge*> narcWidgetEdgeMap;

	// input 
	int nnodes;
	int narcs;
	double* supply;
	int* tail;
	int* head;
	double* obj;
	double* ub;
	double* lb;

	// output
	double* x;
	double* dj;
	double* pi;
	double* slack;

	WidgetGraph* widgetNet;
	
	void constructCplex(WidgetGraph* w){
		widgetNet = w;
		widgetNodeNum = widgetNet->widgetNodes.size();
		widgetEdgeNum = widgetNet->widgetEdges.size();

		int srcOutNum = w->superSource->credit_in_widget_nodes.size() 
			+ w->superSource->debt_in_widget_nodes.size();
		int destInNum = w->superSink->credit_out_widget_nodes.size() 
			+ w->superSink->debt_out_widget_nodes.size();

		nnodes = widgetNodeNum + 2;
		narcs = widgetEdgeNum + srcOutNum + destInNum;

		supply = new double [nnodes];
		tail = new int [narcs];
		head = new int [narcs];
		obj = new double [narcs];
		ub = new double [narcs];
		lb = new double [narcs];

		for (int i = 0; i < nnodes; ++i){
			supply[i] = 0;
		}

		int cnt = 0;
		for (auto widgetEdge : widgetNet->widgetEdges){
			head[cnt] = widgetEdge->nodeTo->getGlobalId();
			tail[cnt] = widgetEdge->nodeFrom->getGlobalId();
			obj[cnt] = 0;
			ub[cnt] = widgetEdge->cap;
			lb[cnt] = 0;
			narcWidgetEdgeMap.push_back(widgetEdge);
			cnt++;
		}

		// src
		supply[nnodes - 2] = w->payment;
		for (auto widgetNodePair : w->superSource->credit_in_widget_nodes){
			head[cnt] = widgetNodePair.second->getGlobalId();
			tail[cnt] = nnodes - 2;
			obj[cnt] = 0;
			ub[cnt] = widgetNodePair.second->originAtomicEdge->capacity;
			lb[cnt] = 0;
			cnt++;
		}
		for (auto widgetNodePair : w->superSource->debt_in_widget_nodes){
			head[cnt] = widgetNodePair.second->getGlobalId();
			tail[cnt] = nnodes - 2;
			obj[cnt] = 0;
			ub[cnt] = widgetNodePair.second->originAtomicEdge->capacity;
			lb[cnt] = 0;
			cnt++;
		}

		// dest
		supply[nnodes - 1] = - w->payment;
		for (auto widgetNodePair : w->superSink->credit_out_widget_nodes){
			head[cnt] = nnodes - 1;
			tail[cnt] = widgetNodePair.second->getGlobalId();
			obj[cnt] = 0;
			ub[cnt] = widgetNodePair.second->originAtomicEdge->capacity;
			lb[cnt] = 0;
			cnt++;
		}
		for (auto widgetNodePair : w->superSink->debt_in_widget_nodes){
			head[cnt] = nnodes - 1;
			tail[cnt] = widgetNodePair.second->getGlobalId();
			obj[cnt] = 0;
			ub[cnt] = widgetNodePair.second->originAtomicEdge->capacity;
			lb[cnt] = 0;
			cnt++;
		}


		x = new double[narcs];
		dj = new double[narcs];
		pi = new double[nnodes];
		slack = new double[nnodes];
	}


	void printInput(){
		for (int i = 0; i < nnodes; ++i){
			cout << "node: " << i << " supply: " << supply[i] << endl;
		}

		for (int i = 0; i < narcs; ++i){
			cout << "head: " << head[i] << " tail: " << tail[i] 
			<< " obj: " << obj[i] << " ub: " 
			<< ub[i] << " lb: " << lb[i] << endl;
		}
	}
	void printResult(){
		for (int i = 0; i < narcs; ++i){
			cout << "head: " << head[i] << " tail: " << tail[i]
			<< " flow: " << x[i] << " dj: " << dj[i]
			<< " pi: " << pi[i] << " slack: " << slack[i] << endl;
		}
	}

	void copyBack(){
		
		for (int i = 0; i < narcWidgetEdgeMap.size(); ++i){
			// cout << "copy back: " << head[i] << " " << tail[i] << " " << x[i] << endl;
			narcWidgetEdgeMap[i]->curr = x[i];
		}

		delete [] supply;
		delete [] head;
		delete [] tail;
		delete [] obj;
		delete [] ub;
		delete [] lb;

		delete [] x;
		delete [] dj;
		delete [] pi;
		delete [] slack;
	}
};


class CplexSolver{

public:

	int solve(CplexConverter& converter);

};


#endif