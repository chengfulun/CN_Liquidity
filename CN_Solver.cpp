#include "CN_Solver.h"
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <climits>
#include <cfloat>

extern CredNetConstants credNetConstants;


bool LpSolver::solveLpProblem(CplexConverter& cplexConverter, string mode, string purpose, double cRate = 0.0, double dRate = 0.1, double haircut = 0.8)
{
	// clock_t t;
	// t = clock();
	
	IloEnv   env;
	bool success = false;
	IloModel model(env);
	IloNumVarArray var(env);
	IloRangeArray con(env);
		try {

		this->populatebyrow(cplexConverter, model, var, con, mode, purpose, cRate, dRate, haircut); //modify where the graph is located
		this->addObjective(mode, cplexConverter, model, var, con);
		IloCplex cplex(model);

		// cout << "before solving " << endl;
		// cout << var << endl;
		// cout << con << endl;

		cplex.setOut(env.getNullStream());
		// Optimize the problem and obtain solution.
		if ( !cplex.solve() ) {
		   // env.error() << "Failed to optimize LP" << endl;
		   throw(-1);
		}

		IloNumArray vals(env);
		
		// env.out() << "Solution status = " << cplex.getStatus() << endl;
		// env.out() << "Solution value  = " << cplex.getObjValue() << endl;
		// cplex.getValues(vals, var);
		// env.out() << "Values        = " << vals << endl;
		// cplex.getSlacks(vals, con);
		// env.out() << "Slacks        = " << vals << endl;
		// cplex.getDuals(vals, con);
		// env.out() << "Duals         = " << vals << endl;
		// cplex.getReducedCosts(vals, var);
		// env.out() << "Reduced Costs = " << vals << endl;
		
    	// cplex.exportModel("lpex1.lp");

		for (int i = 0; i < cplexConverter.variables.size(); ++i){
			cplexConverter.results.push_back(cplex.getValue(var[i]));
		}

		// cout << "after solving " << endl;
		// cout << var << endl;
		// cout << con << endl;
		success = true;

	}
	catch (IloException& e) {
		cerr << "Concert exception caught: " << e << endl;
	}
	catch (...) {
		cerr << "Unknown exception caught" << endl;
	}
	// if(not success){
	// 	cout << var << endl;
	// 	cout << con << endl;
	// }
	env.end();


	// t = clock() - t;
	// printf ("It took me %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);

	return success;
}

void LpSolver::addObjective(string mode, 
	CplexConverter& cplexConverter, IloModel model, 
	IloNumVarArray x, IloRangeArray c){
	
	// cerr << mode << endl;

	IloEnv env = model.getEnv();

	if (mode == "MIN_SRC_COST"){

		IloExpr cost(env);

		// add cost to all atomic edges
		for (int i = 0; i < cplexConverter.variables.size(); ++i){
			cost += x[i] * 0.000000001;
			// cplexConverter.variables[i].interest_rate;
		}

		for(auto &atoIn : cplexConverter.src->atomicEdge_in){
			int aeId = atoIn.second->atomicEdgeId;
			for (int j = 0; j < cplexConverter.atomicIdToVarIdDict[aeId].size(); j++){
				// var Id
				int vId = cplexConverter.atomicIdToVarIdDict[aeId][j];
				cost += cplexConverter.variables[vId].value_paying * x[vId];

				// cout << "adding " << cplexConverter.variables[vId].interest_rate 
				// 		<< " * " << vId << endl;
			}
		}
		model.add(IloMinimize(env, cost));
	}
	else if (mode == "MAX_FLOW_SRC"){

		IloExpr cost(env);

		// add cost to all atomic edges
		for (int i = 0; i < cplexConverter.variables.size(); ++i){
			cost += x[i] * 0.000000001;
			// cplexConverter.variables[i].interest_rate;
		}

		for(auto &atoIn : cplexConverter.src->atomicEdge_in){
			int aeId = atoIn.second->atomicEdgeId;
			for (int j = 0; j < cplexConverter.atomicIdToVarIdDict[aeId].size(); j++){
				// var Id
				int vId = cplexConverter.atomicIdToVarIdDict[aeId][j];
				cost += 1000 -x[vId];
				// cout << "adding " << cplexConverter.variables[vId].interest_rate 
				// 		<< " * " << vId << endl;
			}
		}
		for(auto &atoIn : cplexConverter.src->atomicEdge_in){
			int aeId = atoIn.second->atomicEdgeId;
			for (int j = 0; j < cplexConverter.atomicIdToVarIdDict[aeId].size(); j++){
				// var Id
				int vId = cplexConverter.atomicIdToVarIdDict[aeId][j];
				cost += cplexConverter.variables[vId].value_paying * x[vId];

				// cout << "adding " << cplexConverter.variables[vId].interest_rate 
				// 		<< " * " << vId << endl;
			}
		}		
		model.add(IloMinimize(env, cost));
		cost.end();

	}
	else if (mode == "MAX_FLOW"){

		IloExpr cost(env);

		// add cost to all atomic edges
		for (int i = 0; i < cplexConverter.variables.size(); ++i){
			cost += x[i] * 0.000000001;
			// cplexConverter.variables[i].interest_rate;
		}

		for(auto &atoIn : cplexConverter.src->atomicEdge_in){
			int aeId = atoIn.second->atomicEdgeId;
			for (int j = 0; j < cplexConverter.atomicIdToVarIdDict[aeId].size(); j++){
				// var Id
				int vId = cplexConverter.atomicIdToVarIdDict[aeId][j];
				cost += 1000 -x[vId];
				// cout << "adding " << cplexConverter.variables[vId].interest_rate 
				// 		<< " * " << vId << endl;
			}
		}
		model.add(IloMinimize(env, cost));
		cost.end();

	}	else if (mode == "MAX_FLOW_HIGHIR"){

		IloExpr cost(env);

		// add cost to all atomic edges
		for (int i = 0; i < cplexConverter.variables.size(); ++i){
			cost += x[i] * 0.000000001;
			// cplexConverter.variables[i].interest_rate;
		}

		for(auto &atoIn : cplexConverter.src->atomicEdge_in){
			int aeId = atoIn.second->atomicEdgeId;
			for (int j = 0; j < cplexConverter.atomicIdToVarIdDict[aeId].size(); j++){
				// var Id
				int vId = cplexConverter.atomicIdToVarIdDict[aeId][j];
				cost += 1000 -x[vId]*(1-cplexConverter.variables[vId].interest_rate);
				// cout << "adding " << cplexConverter.variables[vId].interest_rate 
				// 		<< " * " << vId << endl;
			}
		}
		model.add(IloMinimize(env, cost));

	} else if (mode == "MIN_CREDIT_COST") {

		IloExpr cost(env);

		// add cost to all atomic edges
		for (int i = 0; i < cplexConverter.variables.size(); ++i){
			cost += x[i] * cplexConverter.variables[i].interest_rate;
		}

		for (int i = 0; i < cplexConverter.variables.size(); ++i){
			int aeId = cplexConverter.variables[i].atomicEdgeId;
			if (!cplexConverter.graph->atomicEdges[aeId]->isDebt){
				cost += x[i];
			}
		}
		model.add(IloMinimize(env, cost));

	} else if (mode == "MIN_CREDIT_SRC") {

		IloExpr cost(env);

		for(auto &atoIn : cplexConverter.src->atomicEdge_in){
			int aeId = atoIn.second->atomicEdgeId;
			for (int j = 0; j < cplexConverter.atomicIdToVarIdDict[aeId].size(); j++){
				// var Id
				int vId = cplexConverter.atomicIdToVarIdDict[aeId][j];
				cost += cplexConverter.variables[vId].interest_rate * x[vId];
				if (!cplexConverter.graph->atomicEdges[aeId]->isDebt){
					cost += 10*x[vId];
				}

				// cout << "adding " << cplexConverter.variables[vId].interest_rate 
				// 		<< " * " << vId << endl;
			}
		}
		model.add(IloMinimize(env, cost));

	}	else if (mode == "MIN_DEGREE_SRC") {

		IloExpr cost(env);

		for(auto &atoIn : cplexConverter.src->atomicEdge_in){
			int aeId = atoIn.second->atomicEdgeId;
			AtomicEdge* atEdge = cplexConverter.graph->atomicEdges[aeId];			
			for (int j = 0; j < cplexConverter.atomicIdToVarIdDict[aeId].size(); j++){
				// var Id
				int vId = cplexConverter.atomicIdToVarIdDict[aeId][j];
				cost += cplexConverter.variables[vId].interest_rate * x[vId];					
				if (!atEdge->isDebt){
					cost += atEdge->nodeTo->degree * x[vId];
				} else {
					cost += atEdge->nodeFrom->degree * x[vId];
				}
				// cout << "adding " << cplexConverter.variables[vId].interest_rate 
				// 		<< " * " << vId << endl;
			}
		}
		model.add(IloMinimize(env, cost));
		cost.end();

	}	
		else if (mode == "MAX_IR_COST") {

		IloExpr cost(env);
		// add cost to all atomic edges
		for (int i = 0; i < cplexConverter.variables.size(); ++i){
			cost += 1000 - x[i] * cplexConverter.variables[i].interest_rate;
		}
		model.add(IloMinimize(env, cost));

	}	else if (mode == "MIN_SUMIR_COST") {

		IloExpr cost(env);
		// add cost to all atomic edges
		for (int i = 0; i < cplexConverter.variables.size(); ++i){
			cost += x[i] * cplexConverter.variables[i].interest_rate;
		}
		model.add(IloMinimize(env, cost));

	} 	else if (mode == "MIN_DEGREE_COST") {

		IloExpr cost(env);

		// add cost to all atomic edges
		for (int i = 0; i < cplexConverter.variables.size(); ++i){
			cost += x[i] * cplexConverter.variables[i].interest_rate;
		}
		
		for (int i = 0; i < cplexConverter.variables.size(); ++i){
			int aeId = cplexConverter.variables[i].atomicEdgeId;
			AtomicEdge* atEdge = cplexConverter.graph->atomicEdges[aeId];
			
			for (int j = 0; j < cplexConverter.atomicIdToVarIdDict[aeId].size(); j++){
				// var Id
				int vId = cplexConverter.atomicIdToVarIdDict[aeId][j];
				if (!atEdge->isDebt){
					cost += atEdge->nodeTo->degree * x[vId];
				} else {
					cost += atEdge->nodeFrom->degree * x[vId];
				}	
			}

		}
		model.add(IloMinimize(env, cost));

	} else {
		// default
		model.add(IloMinimize(env, 1));
	}


}
	
void LpSolver::populatebyrow (CplexConverter& cplexConverter, 
	IloModel model, IloNumVarArray x, IloRangeArray c, string mode, string purpose, double cRate, double dRate, double haircut)
{
	IloEnv env = model.getEnv();
	// CAPITAL LETTERS MEAN I NEED YOUR HELP, here is help 

	// IloExpr cost(env);
	
	// Create Variables
	// cout << "size of var: " << cplexConverter.variables.size() << endl;
	for (int i = 0; i < cplexConverter.variables.size(); ++i){
		IloNumVar iloVar(env, 0.0, cplexConverter.capacities[i]);
		// cout << iloVar << endl;
		x.add(iloVar);
	}

	// Capacity Constraints
	for (auto &it : cplexConverter.atomicIdToVarIdDict){
		IloExpr t(env);
		// cout << "adding constraint ";
		// loop below for credit edges that have multiple variables
		for (int j = 0; j < it.second.size(); j++){
			// cout << "x[" << it.second[j] << "] + ";
			t += x[it.second[j]];
		}
		// cout << endl;
		c.add(t <= cplexConverter.graph->atomicEdges[it.first]->capacity);
		// cout << c << endl;
		t.end();
	}

	// other constraints
	for (auto nodePair : cplexConverter.graph->nodes){

		// For all nodes
		Node* n = nodePair.second;

		if(n == cplexConverter.src){

			// source constraints
			// IloExpr inFlow(env);
			IloExpr outFlow(env);
			IloExpr collateral_req(env);
			IloExpr collateral(env);

			for(auto &atoIn : n->atomicEdge_in){
				int aeId = atoIn.second->atomicEdgeId;
				int ceId = atoIn.second->singleCreditIndex;
				double cr = atoIn.second->originEdge->singleCreditEdges[ceId]->collateralRate;
				// cout<<"cr "<<cr<<endl;			
				for (int j = 0; j < cplexConverter.atomicIdToVarIdDict[aeId].size(); j++){
					// var Id
					int vId = cplexConverter.atomicIdToVarIdDict[aeId][j];
					outFlow += x[vId];
					if (not atoIn.second->isDebt){
						collateral_req += x[vId] * cr;
					}

					if(purpose == "DEBT"){
						// input cRate at the time of calling payAsset()
						collateral_req -= x[vId] * cRate;
					}

					if(atoIn.second->nodeTo->isMarket && atoIn.second->isDebt && purpose == "PURCHASE"){
						collateral -= x[vId];
					}
					if(purpose == "ASSET"){
						collateral -= x[vId] * (1-haircut);
						// collateral -= x[vId] * haircut;
					}

					// cost += cplexConverter.graph->atomicEdges[cplexConverter.variables[vId].atomicEdgeId]->interest_rate * x[vId];
				}
			}

			if (not n->isMarket){
				// set dRate in payAsset
				collateral_req += n->getCollateral(dRate);
				collateral += n->getWealth(haircut);
				// cout<<"wealth added "<<n->getWealth(haircut);
			}
			// make sure no flow goes into source
			for (auto &atoOut : n->atomicEdge_out){
				int aeId = atoOut.second->atomicEdgeId;

				for (int j = 0; j < cplexConverter.atomicIdToVarIdDict[aeId].size(); j++){
					// var Id
					int vId = cplexConverter.atomicIdToVarIdDict[aeId][j];
					// inFlow += x[vId];
					c.add(x[vId] == 0);
					// cost -= cplexConverter.graph->atomicEdges[cplexConverter.variables[vId].atomicEdgeId]->interest_rate * x[vId];
				}
			}
			if (mode == "MAX_FLOW" || mode == "MAX_FLOW_SRC" || purpose == "ASSET"){
				c.add(outFlow <= cplexConverter.request);
			}
			else{
				// c.add(outFlow == cplexConverter.request);				
			}

			if(not n->isMarket and not n->leveraged && not n->defaulted && n->getWealth(haircut) > n->getCollateral(dRate)){
				c.add(collateral  - collateral_req >= 0);
			}
			// inFlow.end();
			collateral.end();
			collateral_req.end();
			outFlow.end();

		} else if(n == cplexConverter.dest){

			// destination constraints
			IloExpr inFlow(env);
			// IloExpr outFlow(env);
			// make sure dest doesn't pay anything
			for(auto &atoIn : n->atomicEdge_in){
				int aeId = atoIn.second->atomicEdgeId;
				for (int j = 0; j < cplexConverter.atomicIdToVarIdDict[aeId].size(); j++){
					// var Id
					int vId = cplexConverter.atomicIdToVarIdDict[aeId][j];
					// outFlow += x[vId];
					c.add(x[vId] == 0);
				}
			}
			for (auto &atoOut : n->atomicEdge_out){
				int aeId = atoOut.second->atomicEdgeId;
				for (int j = 0; j < cplexConverter.atomicIdToVarIdDict[aeId].size(); j++){
					// var Id
					int vId = cplexConverter.atomicIdToVarIdDict[aeId][j];
					inFlow += x[vId];
				}
			}
			if (mode == "MAX_FLOW" || mode == "MAX_FLOW_SRC" || purpose == "ASSET"){
				c.add(inFlow <= cplexConverter.request);
			}
			else{
				// c.add(inFlow == cplexConverter.request);				
			}
			inFlow.end();
			// outFlow.end();

		} else {

			// Monotonicity Constraints
			for (int i = 0; i < credNetConstants.totalValues.size(); ++i){
				IloExpr tempin(env);
				IloExpr tempout(env);

				for (auto &atoIn : n->atomicEdge_in){
					int aeId = atoIn.second->atomicEdgeId;
					for (int j = 0; j < cplexConverter.atomicIdToVarIdDict[aeId].size(); j++){

						// var Id
						int vId = cplexConverter.atomicIdToVarIdDict[aeId][j];
						// if(i==0 && cplexConverter.variables[vId].value_paying==credNetConstants.totalValues[0]){
						// 	tempout += x[vId];
						// }							
						if (cplexConverter.variables[vId].value_paying <= credNetConstants.totalValues[i]){
							tempout += x[vId];
						}
					}
				}
				for (auto &atoOut : n->atomicEdge_out){
					int aeId = atoOut.second->atomicEdgeId;
					for (int j = 0; j < cplexConverter.atomicIdToVarIdDict[aeId].size(); j++){

						// var Id
						int vId = cplexConverter.atomicIdToVarIdDict[aeId][j];
						// if(i==0 && cplexConverter.variables[vId].interest_rate == credNetConstants.totalIrs[0]){
						// 	tempin += x[vId];							
						// }
						if (cplexConverter.variables[vId].value_receiving <= credNetConstants.totalValues[i]){
							tempin += x[vId];
						}
					}
				}
				// if(i==0){
				// 	// c.add(tempin <= DBL_MIN);					
				// }
				// if(i>2){
				if(not n->isMarket){
					c.add(tempout + tempin >= 0);						
				}
				// }
				tempout.end();
				tempin.end();
			}

			//Flow Constraints
			IloExpr inFlow(env);
			IloExpr outFlow(env);
			IloExpr collateral_req(env);
			IloExpr collateral(env);
			for(auto &atoIn : n->atomicEdge_in){
				int aeId = atoIn.second->atomicEdgeId;
				int ceId = atoIn.second->singleCreditIndex;
				double cr = atoIn.second->CR();
				// originEdge->singleCreditEdges[ceId]->collateralRate;
				for (int j = 0; j < cplexConverter.atomicIdToVarIdDict[aeId].size(); j++){
					// var Id
					int vId = cplexConverter.atomicIdToVarIdDict[aeId][j];
					outFlow += x[vId];
					if (not atoIn.second->isDebt){
						collateral_req += x[vId] * cr;
					}
					// if (atoIn.second->nodeTo->isMarket && atoIn.second->isDebt && purpose=="PURCHASE"){
					// 	collateral -= x[vId];
					// }

				}
			}
			if (not n->isMarket){
				collateral_req += n->getCollateral(dRate);
			}
			for (auto &atoOut : n->atomicEdge_out){
				int aeId = atoOut.second->atomicEdgeId;
				int ceId = atoOut.second->singleCreditIndex;
				double cr = atoOut.second->CR();
				// originEdge->singleCreditEdges[ceId]->collateralRate;				
				for (int j = 0; j < cplexConverter.atomicIdToVarIdDict[aeId].size(); j++){
					// var Id
					int vId = cplexConverter.atomicIdToVarIdDict[aeId][j];
					inFlow += x[vId];

					// if (atoOut.second->nodeTo->isMarket && not atoOut.second->isDebt){
					// 	collateral += x[vId];
					// }
					if (atoOut.second->isDebt){
						collateral_req -= x[vId] * cr;
					}
				}
			}
			// double paySum = 0;
			// for (auto &toPay : n->creditPayOut){
			// 	paySum += toPay.second.second;
			// }
			if (not n->isMarket){
				collateral += n->getWealth(haircut);
				// collateral -= paySum;			
			}
			// collateral += n->assets;

			c.add(inFlow - outFlow == 0);
			// && not n->defaulted && n->getWealth() > n -> getCollateral()			
			if(not n->isMarket and not n->leveraged && not n->defaulted && n->getWealth(haircut) > n -> getCollateral(dRate)){
				c.add(collateral - collateral_req >= 0);
			}
			inFlow.end();
			outFlow.end();
			collateral.end();
			collateral_req.end();

		}

	}


	model.add(c);
	// model.add(IloMinimize(env, cost));
	// model.add(IloMaximize(env,cost));  //option to minimize cost
	// cost.end();

}  // END populatebyrow

