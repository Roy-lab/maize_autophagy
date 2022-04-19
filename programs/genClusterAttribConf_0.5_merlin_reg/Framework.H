#ifndef _FRAMEWORK_
#define _FRAMEWORK_
#include <map>
#include <string>
using namespace std;
class GeneNameMapper;
class GeneExpManager;
class Framework
{
	public:
	Framework();
	~Framework();
	int init(int, char**);
	int genClusterAttribs();
	int genClusterRegulators();
	private:
	int readSelectClusterEnrPairs(const char*);
	int readEnr(const char*,string&);
	int readClusterMembership(const char*);
	int readHeaderNames(const char*);
	
	// DC added
	int readEdgeConfs(const char*);

	private:
	int minT;
	int maxT;
	int readConfig(const char*);
	int orderColumns(int,map<string,int>*,int,vector<string>&);
	int orderColumns_Optimal(int,map<string,int>*,int,vector<string>&);
	int orderGenes(int,map<string,int>*,int,vector<string>&);
	int orderGenes_Cnt(int,map<string,int>*,int,vector<string>&);
	int orderGenes_Optimal(int,map<string,int>*,int,vector<string>&);
	int orderModAvg_Optimal(map<string, vector<double>* >& modAvgMap, vector<string>& modOrder);
	map<int,map<string,int>*> clusterEnrPair;
	map<string,map<string,string>*> geneEnrMap;
	map<int,map<string,int>*> exprClust;
	GeneNameMapper gnm;
	GeneExpManager expMgr;
	map<int,string> headerNames;
	char outDirName[1024];

	// DC added - conf values for expression regs to targets
	// Key is "REG-TARGET" string; value is read in from file
	map<string,double> edgeConfs;
	bool hasReadEdgeConfs;
	vector<double>* getAverage(vector< vector<double>* >& moduleExp);
};
#endif
