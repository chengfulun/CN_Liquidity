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

	CreditNet(double dVol, int finNumT, double precision, int marketId, 
		double initR,double initVol,double drate, double deposit,
		double haircut, double mReserve, double mLimit,
		int valueBins, double EAR, double asset_vol, int defaulted_periods, double explore_boost, bool verb);
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
	void activateReserves(int fid);
	//mode 0 means only reserve for debt repayment, not lending
	void deactivateReserves(int fid, int mode =0);
	void print_values();
	void init_lambdas();
	// vector<double> returns;
	// vector<double> volatilities;
	// vector<double> wealths;
	// vector<double> credits;
	// vector<double> credits_last;
	vector<double> cReturns;
	vector<double> c2Returns;

	double expected_deposits;
	// credit market
	vector<double> credit_requests;
	vector<double> credit_premiums;
	vector<double> credit_utils;
	vector<double> active_num;
	vector<double> inactive_num;

	vector<double> default_rates;
	// total cost of debt used to fund debt investments
	vector<double> debt_coc;
	// total cost of debt used to pay for assets
	vector<double> asset_coc;
	vector<double> asset_amt;
	vector<double> debt_amt;

	vector< vector<int> > active_set_statuses;

	// initial options
	bool alternate = 0;
	double precision;
	int transactionCounter;
	double mReserve;	

	// for updateReturns()
	double initR;	
	double initVol;
	double total_debt;
	double total_assets;

	double asset_coc_all;
	double debt_coc_all;
	double creturn_all;
	double creturn2_all;

	double mLimit;
	// double dr_all;
	// double cr_all;
	// double ir_all;
	// double ir2_all;
	// double one_minus_cr2_all;
	// double one_minus_dr2_all;
	// double dr2_all;

	vector<int> rDefaults;
	vector<int> cDefaults;
	vector<int> dDefaults;
	vector<int> aDefaults;
	vector<int> fDefaults;
	vector<int> liquid_loss;
	vector<int> default_wealth;
	vector<int> deposit_losses;
	vector<int> creditor_losses;
	vector<double> all_values;

	void record_values();

	int shockPay(double alpha, bool verbose);
	double getReturn(int nodeNum, string mode);
	int ntrials = 0;
	void updateReturns();
	int makeInvest(bool v, bool verb);
	double makeLiquidate(int fid, double amt);
	vector<double> credit_shares(int ix);

	std::tuple<double,double,double> payAsset(int fid1, int fid2, double amt, string mode, int transSeqNum, string purpose, double drate, double haircut, double crate);
	// (check,asset_coc,debt_coc)
	void updateValueBins();
	double deposit_rate;
	double haircut;
	// vector<double> credits_sq;
	// vector<double> cReturns_sq;
	// vector<double> credits_cross;
	int maturity;
	void resultsOut();
	void resultsOut_1();
	bool verb;
	int fDefault= 10;
	int cCount= 12;
	vector<double> defaultRates;
	// double liquid_loss = 0.0;
	// double default_wealth = 0.0;
	// double deposit_losses =0.0;
	// double creditor_losses = 0.0;
	vector<double> assetPrices;
	int valueBins;

	double price;
	double price_last;
	void setReserves();
	double collateralValue(int nodeIx, double CR, double price, bool useCredit);
	void updateCollateralValues();
	void requestLines();
	void update_default_rates();
	double defaultRate(int nodeIx);
	void updateLineRequests();
	void updateCreditUtil();
	void updateLambdas();
	void updateCreditPremiums();
	void updateCreditTerms();
	void updateCreditLimits();
	void updateCreditSets(double thresh, int inactive_periods);
	void updateCreditPremiums(double thresh, double uprate, double downrate);
	void restore_defaulted(int fid);
	void postCashUpdateLambda();

	vector<double> debt2_amt;

	void makeDefault(int fid);
	void clear_history();
	double expected_asset_return;
	int defaulted_periods;
	double asset_volatility;
	double explore_boost;

	double log_mu;
	double log_vol;
	double debt_coc_total;
	vector<double> cReturns_real;
	vector<double> debt_vol;
	double alpha;
	double payCash(int fid, double amt);
	double dVol;
	double dRevert;

};


#endif