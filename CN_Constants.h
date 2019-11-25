#ifndef CN_Const
#define CN_Const

#include <vector>
#include <random>
#include <algorithm>
#include <chrono>
#include <set>

using namespace std;

class CredNetConstants{

public:
	vector<double> totalIrs;
	double assetPrice;
	vector<double> totalCrs;
	vector<double> totalDrs;
	vector<double> totalArs;
	vector<double> totalValues;

	default_random_engine globalGenerator;
	uniform_int_distribution<int> uniformIntDistribution;
	uniform_real_distribution<double> uniformDoubleDistribution;
	normal_distribution<double> normalDistribution;


	CredNetConstants(): uniformIntDistribution(0, 9999), uniformDoubleDistribution(0.0, 1.0)
		, normalDistribution(0.0,1.0),globalGenerator(std::chrono::system_clock::now().time_since_epoch().count()){
	}

	void addIr(double ir);
	void addCr(double cr);
	void addDr(double dr);
	void addAr(double ar);
	void setAllValues(vector<double> irs,vector<double> crs,vector<double> drs,vector<double> ars);
	void clean();
	void priceUpdate(double price);
	void print();
	void setValues(vector<double> values, int k);
};


#endif