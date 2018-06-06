	
#include "CN_Constants.h"
#include <iostream>
#include <random>

using namespace std;

CredNetConstants credNetConstants;

void CredNetConstants::print(){
	cout << "possible irs: ";
	for (int i = 0; i < totalIrs.size(); ++i){
		cout << totalIrs[i]<<" ,";
	}
	cout << "possible crs: ";
	for (int i = 0; i < totalCrs.size(); ++i){
		cout << totalCrs[i]<<" ,";
	}
	cout << endl;
}


void CredNetConstants::addIr(double ir){
	totalIrs.push_back(ir);
}

void CredNetConstants::addCr(double cr){
	totalCrs.push_back(cr);
}

void CredNetConstants::priceUpdate(double val){
	assetPrice = val;
}

void CredNetConstants::clean(){
	totalIrs.clear();
    totalCrs.clear();
}