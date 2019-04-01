#include "include/rapidjson/document.h"     // rapidjson's DOM-style API
#include "include/rapidjson/prettywriter.h" // for stringify JSON
#include "include/rapidjson/filewritestream.h"
#include "include/rapidjson/filereadstream.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <thread>
#include <algorithm>
#include <mutex>
#include <random>
#include <iomanip>

using namespace rapidjson;
using namespace std;

#include "CN_CreditNet.h"
extern CredNetConstants credNetConstants;
class Generator {
	std::default_random_engine generator;
	std::normal_distribution<double> distribution;
	double min;
	double max;
public:
	Generator(double mean, double stddev, double min, double max):
		distribution(mean, stddev), min(min), max(max)
	{}

	double operator ()() {
		while (true) {
			double number = this->distribution(generator);
			if (number >= this->min && number <= this->max)
				return number;
		}
	}
};

struct Config {
	string numNodes;
	string smoothing;
	string numIR;
	string epochs;
	string FFR;
	string wealth;
	string deposits;
	string dShock;
	string CR;
	string EAR;
	string asset_vol;
	string dRate;
	string haircut;
	string initR;
	string initV;
	string maturity;
	vector<string> assignedStrategy;
};

struct PlayerInfo {
	double payoff;
	string strategy;
	string role;
	// int features[10];
};

void readConfig (Config &config, string inPath) {
	// cout<<"trying to read"<<endl;
	FILE* fp = fopen(inPath.c_str(), "r"); // non-Windows use "r"
	// cout<<"fopened"<<endl;
	char readBuffer[65536];
	// cout << "buffer set" << endl;
	FileReadStream is(fp, readBuffer, sizeof(readBuffer));
	// cout << "stream set" << endl;
	Document doc;
	doc.ParseStream(is);
	fclose(fp);
	// cout<<"fclosed"<<endl;

	// printf("\nModified JSON with reformatting:\n");
	// StringBuffer sb;
	// PrettyWriter<StringBuffer> writer(sb);
	// doc.Accept(writer);    // Accept() traverses the DOM and generates Handler events.
	// puts(sb.GetString());
	
	const Value& configObj = doc["configuration"];
	config.smoothing = configObj["smoothing"].GetString();
	config.numNodes = configObj["numNodes"].GetString();
	config.numIR = configObj["numIR"].GetString();
	config.epochs = configObj["epochs"].GetString();
	config.FFR = configObj["FFR"].GetString();
	config.wealth = configObj["wealth"].GetString();
	config.deposits = configObj["deposits"].GetString();
	config.dShock = configObj["dShock"].GetString();
	config.CR = configObj["CR"].GetString();
	config.EAR = configObj["EAR"].GetString();
	config.asset_vol = configObj["asset_vol"].GetString();
	config.dRate = configObj["dRate"].GetString();
	config.haircut = configObj["haircut"].GetString();
	config.initR = configObj["initR"].GetString();
	config.initV = configObj["initV"].GetString();
	config.maturity = configObj["maturity"].GetString();
	

	const Value& a = doc["assignment"];
	const Value& b = a["All"];
	
	// rapidjson uses SizeType instead of size_t.
	for (rapidjson::SizeType i = 0; i < b.Size(); i++)
	{
		config.assignedStrategy.push_back(b[i].GetString());
		// printf("%s \n", b[i].GetString());
		// cout << b[i].GetString() << endl;
	}
}

void writePayoff (std::vector<PlayerInfo> &players, string outPath) {
	
	rapidjson::Document result;
	result.SetObject();
	rapidjson::Value playerArray(rapidjson::kArrayType);
	rapidjson::Document::AllocatorType& allocator = result.GetAllocator();
	
	for (int i = 0; i < players.size(); ++i) {
		// create a rapidjson object type
		rapidjson::Value s;
		s = rapidjson::StringRef(players[i].strategy.c_str());
		rapidjson::Value object(rapidjson::kObjectType);
		object.SetObject();
		object.AddMember("role", "All", allocator);
		// cout<<players[i].strategy<<endl;
		// cout<<players[i].strategy.c_str()<<endl;
		object.AddMember("strategy", s, allocator);
		object.AddMember("payoff", players[i].payoff, allocator);
		playerArray.PushBack(object, allocator);
	}
	
	result.AddMember("players", playerArray, allocator);
	
	rapidjson::Value object(rapidjson::kObjectType);
	object.SetObject();
	result.AddMember("features", object, allocator);
	// printf("\nModified JSON with reformatting:\n");
	// StringBuffer sb;
	// PrettyWriter<StringBuffer> writer(sb);
	// result.Accept(writer);    // Accept() traverses the DOM and generates Handler events.
	// puts(sb.GetString());
	
	FILE* fp = fopen(outPath.c_str(), "w"); // non-Windows use "w"
	char writeBuffer[65536];
	FileWriteStream os(fp, writeBuffer, sizeof(writeBuffer));
	Writer<FileWriteStream> writer1(os);
	result.Accept(writer1);
	fclose(fp);
}

int main(int argc, char* argv[]){
	string json_folder = argv[1];
	// cout<<json_folder<<endl;
	int num_obs = atoi(argv[2]);
	Config config;


	// move to config file in num_obs loop
	// double maxCR = atof(argv[3]);
	// int epochs = atoi(argv[4]);
	// Generator g(2.0, 1.0, 0.0, 10.0);
	// cout <<  num_obs << endl;
	double precision = 100000;
	readConfig(config, json_folder+"/csimspecFDH.json");
			// cout << "configed" << endl;

	int numIR = atoi(argv[3]);
	// cout<<"read config"<<endl;
 	double FFR = stod(config.FFR);

	credNetConstants.clean();

	for (int i = 0; i < numIR; i++){
		// credNetConstants.addIr(i);
		if(i==0){
			credNetConstants.addIr(0.0);
		}
		else{
			credNetConstants.addIr(FFR + (i-1)/200.0);			
		}
	}

	// double threshold = stod(config.edgeProb);
	// double threshold = atof(argv[5]);
	int iter = stoi(config.smoothing);
	int finNum = stoi(config.numNodes);
	int epochs = stoi(config.epochs);
	std::vector<double> payoffs(finNum - 1,0.0);
 	// vector<double> capacities = {10};
 	// {3,5,7,9,11}
 	// *precision;
 	// *precision;
 	// double wealth = atof(argv[4]);
 	double wealth = stod(config.wealth);
 	double deposit = stod(config.deposits);
 	double shock = stod(config.dShock);
 	double cr = stod(config.CR);
 	shock = atof(argv[4]); 	
 	double EAR = stod(config.EAR);
 	double asset_vol = sqrt(stod(config.asset_vol));
 	double deposit_rate = stod(config.dRate);
 	double haircut = stod(config.haircut);
 	double initR = stod(config.initR);
 	double initV = stod(config.initV); 	
 	int maturity = stoi(config.maturity); // minus one
 	bool outV = false;
 	bool outResults = true;
 	bool randThetas = false;
 	int forceT = epochs - 5;
 	// {1,2,3,4,5}
 	// * precision;
 	// cout<<"capacity is "<<capacity<<endl;
 	// cout<<"transaction amount is "<<transAmount<<endl;
			for (int i = 0; i < num_obs; ++i){
				// cout<<"ob started"<<endl;
				for (int j = 0; j < iter; ++j){
					int failRateTotal = 0;
					int rDefaults = 0;
					int cDefaults = 0;
					int rDefaultsAlt;
					int cDefaultsAlt;					
					int marketId = finNum - 1;
					CreditNet creditNet(finNum,precision, marketId, initR, initV, deposit_rate, haircut, maturity, outV);
					// cout<<"initialized crednet"<<endl;
					creditNet.genMarket0Graph(deposit, shock, wealth, FFR, marketId, cr, EAR, asset_vol, randThetas);
					// creditNet.genTest0Graph(threshold, numIR, capacity,maxCR,wealth,marketId);
								// cout<<"generated graph"<<endl;
					// creditNet.print();

					creditNet.setThetas(config.assignedStrategy);
					// CreditNet* creditNetAlt = new CreditNet();
					CreditNet creditNetAlt(finNum,precision, marketId, initR, initV, deposit_rate, haircut, maturity, outV);

					creditNet.alternate = false;
					for (int i_e = 0; i_e < epochs; ++i_e){
						bool vv = false;
						bool fd = false;
						if(i_e == forceT){
							fd = true;
							creditNetAlt = creditNet;
							creditNetAlt.alternate = true;
							rDefaultsAlt = rDefaults;
							cDefaultsAlt = cDefaults;					
						}
						// if((i_e+1)%10==0 && i_e>0 && outResults){
						// 	vv = true;
						// }				
						if((i_e >forceT - 2) && i_e>0 && outResults){
							vv = true;
						}							
						// double price = g();
						// cout<<"start"<<endl;
						if(outV){
							cout<<"Epoch----------------"<<i_e<<endl;
						}
						cDefaults += creditNet.makeInvest(fd,vv);
						// creditNet.print();
						double rr = (credNetConstants.normalDistribution(
											credNetConstants.globalGenerator)*asset_vol*sqrt((double)(maturity - 1)));
						double alpha =  exp((EAR - asset_vol*asset_vol/2.0)*(double)(maturity - 1)+rr);
						if(outV){cout<<"shock done "<<alpha<<endl;}
						// cout << window_size - failRateTotal << "   "<<endl;
						rDefaults += creditNet.shockPay(alpha, false);

						if(i_e>=forceT){
							cDefaultsAlt += creditNetAlt.makeInvest(false,vv);
							rDefaultsAlt += creditNetAlt.shockPay(alpha, false);
							creditNetAlt.fDefaults = creditNet.fDefaults;
						}
							// cout<<"paid IR"<<endl;
							// cout<<"node "<< l << " processed"<<endl;
							// cDefaults += creditNet.checkCollateral(l);
							// cout<<"paid cr"<<endl;
						// creditNet.print();
						// cout<<"epoch "<<i_e<<" done"<<endl;
							 // cout  << payoffs[k] << "   ";

						// creditNet.nodes[finNum-1]->print();
						// cout << numIR << "  "<< deposit  <<"  "<<rDefaults  <<endl;
						// cout << "failed trx: " << failRateTotal << endl;
						// cout << "rDefaults: " << rDefaults << endl;
						// cout << "cDefaults: " << cDefaults << endl;
						// d += creditNet.checkCR(price);
						// defaults += d;
					
						// cout<<"payoffs written"<<endl;
						// creditNet.print();
												
					

				// cout<<"exit loop"<<endl;
					
				
			}
				for (int k = 0; k < finNum - 2; k++){
							// cout << creditNet.nodes[k]->transactionNum << "  " << creditNet.nodes[k]->getCurrBalance()/(precision*100)<<"   ";
							payoffs[k] += (double) creditNet.nodes[k]->getWealth(1.0);
								// *precision
							// cout << endl;
						}

				// for (int k = 0; k < finNum - 1; ++k){
				// 			// cout << creditNet.nodes[k]->transactionNum << "  " << creditNet.nodes[k]->getCurrBalance()/(precision*100)<<"   ";
				// 	delete creditNet.nodes[k];								// *precision
				// 			// cout << endl;
				// 		}						

		}
	// cout<<"exit loop 2"<<endl;
						// cout<<"payouts ready"<<endl;
						std::vector<PlayerInfo> myList;
						for (int g = 0; g < finNum - 2; ++g) {
							PlayerInfo p;
							p.strategy = config.assignedStrategy[g];
							p.payoff = (double)payoffs[g]/(double)iter;
							p.role = "All";
							myList.push_back(p);
							}
						// cout<<"payoffs ready"<<endl;

						writePayoff(myList, json_folder + "/observation" + std::to_string(i) + ".json");		
	}

	return 0;	
}