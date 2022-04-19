#ifndef _SUBGRAPH_MGR
#define _SUBGRAPH_MGR

#include <vector>
using namespace std;

class SubGraph;

class SubGraphMgr
{
	public:
		SubGraphMgr();
		~SubGraphMgr();
		int readSubGraphs(const char*);
		int readSubGraphAttributes(const char*);
		int readGeneList(const char*);
		int genRandomGraphs(int,int);
		int genRandomGraphs(int,int,vector<SubGraph*>&);
		vector<SubGraph*>& getSubGraphs();
		map<string,int>& getAttributeNames();
		SubGraph* getGraphAt(int);
		int getGeneCnt();
		map<string,int>& getGeneList();
		string& getSGName(int);
	private:
		vector<SubGraph*> subGraphSet;
		map<string,int> geneNameIDMap;
		map<int,string> geneIDNameMap;
		map<string,int> attributeNames;
		map<int,string> idNameMap;

};
#endif
