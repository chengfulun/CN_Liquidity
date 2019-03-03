#ifndef CN_CreditNet
#define CN_CreditNet

#include <unordered_map>
#include <list>
#include <string>
#include "CN_Node.h"
#include "CN_Graph.h"
#include "CN_Solver.h"

class CreditNet : public Graph{
private:

public:

	bool alternate = 0;
	double precision;
	int transactionCounter;
	double initVol;
	CreditNet(int finNumT, double precision, int marketId, double initR,double initVol,double drate, double haircut, int maturity, bool verb);
	CreditNet();
	CreditNet(CreditNet &graphT);
	CreditNet& operator=(CreditNet &graphT);
	// Inter Bank Trans
	// liquidity test
	int genInterBankTrans(double request, string mode, int transSeqNum);
	int payIR(int fid1);
	int checkCollateral(int fid1);
	void printPayoff();
	void unwind(int fid);
	void grabCollateral(int fid, double amt);
	int pay(int fid1, int fid2, double amt, string mode, int transSeqNum);
	int payCollateral(int fid, double amt);
	vector<double> returns; // I hope this is CR in 3.6
	vector<double> volatilities; // I hope this is CV in 3.6
	vector<double> wealths;
	vector<double> credits;
	vector<double> credits_last; // how much credit at the beginning of the period
	vector<double> cReturns; // possibly historical return 
	vector<int> rDefaults;
	vector<int> cDefaults;
	vector<int> dDefaults;
	vector<int> aDefaults;
	vector<int> fDefaults;
	vector<int> liquid_loss;
	vector<int> default_wealth;
	vector<int> deposit_losses;
	vector<int> creditor_losses;
	int shockPay(double alpha, bool verbose);
	int ntrials = 0;
	void updateReturns();
	int makeInvest(bool v, bool verb);
	double makeLiquidate(int fid, double amt);
	vector<double> credit_shares(int ix);
	double payAsset(int fid1, int fid2, double amt, string mode, int transSeqNum, string purpose, double crate, double drate, double haircut);
	double deposit_rate;
	double haircut;
	// vector<double> credits_sq;
	// vector<double> cReturns_sq;
	// vector<double> credits_cross;
	int maturity;
	double initR;
	void resultsOut();
	void resultsOut_1();
	bool verb;
	int fDefault= 10;
	int cCount= 12;
	// double liquid_loss = 0.0;
	// double default_wealth = 0.0;
	// double deposit_losses =0.0;
	// double creditor_losses = 0.0;
};


#endif