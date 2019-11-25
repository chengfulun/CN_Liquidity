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

#include "CN_Constants.h"


int main(int argc, char* argv[]){

	CredNetConstants credNetConstants;
	credNetConstants.clean();
	vector<double> testV{1,2,34,34,34,34,87,88};
	cout<<"input"<<endl;
	for (auto &v : testV){
		cout<<v<<endl;
	}
	credNetConstants.setValues(testV, 3);

	return 0;	
}