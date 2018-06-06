#include "CN_Executer.h"

#include <vector>
#include <algorithm>

void Executer::execute(vector<Edge*>& path, Node* src, Node* dest){
	if (path.empty()){
		return;
	}

	vector<bool> isDebt;
	for (int i = 0; i < path.size(); ++i){
		isDebt.push_back(path[i]->get_d_current() == 1);
	}


	int i = path.size() - 1;
	double currIR = isDebt[i] ? 
		path[i]->get_debt_interest_rate() : path[i]->get_interest_rate();
	if (isDebt[i]){
		path[i]->set_debt_interest_rate(path[i]->get_interest_rate());
	}
	path[i]->setCurrent(isDebt[i] ? 0 : 1);
	--i;

	for (; i >= 0; --i){
		path[i]->setCurrent(isDebt[i] ? 0 : 1);

		if (isDebt[i] && isDebt[i+1]){

			currIR = path[i]->get_debt_interest_rate();
			path[i]->set_debt_interest_rate(path[i]->get_interest_rate());

		} else if (!isDebt[i] && isDebt[i+1]) {

			currIR = max(path[i]->get_interest_rate(), currIR);
			path[i]->set_debt_interest_rate(currIR);

		} else if (isDebt[i] && !isDebt[i+1]) {

			currIR = path[i]->get_debt_interest_rate();
			path[i]->set_debt_interest_rate(path[i]->get_interest_rate());

		} else if (!isDebt[i] && !isDebt[i+1]) {

			currIR = max(path[i]->get_interest_rate(), currIR);
			path[i]->set_debt_interest_rate(currIR);

		}

	}

	src -> transactionNum++;
	src -> totalIR += currIR;

}