#ifndef _FRAMEWORK_
#define _FRAMEWORK_
#include <iostream>
#include <map>
#include <string>
#include <vector>
using namespace std;
class HierarchicalCluster;
class OptimalLeafOrder;
class Framework
{
	public:
	Framework();
	~Framework();
	int readDataMatrix(const char*);
	int readDataPair(const char*);
	int reorder(const char*);
	int reorder(const char*,double);
	private:
	int getTree();
	int optimalLeafOrder();
	int readColumns(char*);
	HierarchicalCluster cluster;
	OptimalLeafOrder olo;
	map<string,HierarchicalClusterNode*> nodeSet;
	map<string,HierarchicalClusterNode*> backup;
	map<string,int> nameIDMap;
	map<string,int> attribNameIDMap;
	map<int,string> attribIDNameMap;
};
#endif
