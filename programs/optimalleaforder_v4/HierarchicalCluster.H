#ifndef _HIERARCHICAL_CLUSTER
#define _HIERARCHICAL_CLUSTER
#include <map>
#include <string>
#include <list>
#include <vector>
#include <queue>
using namespace std;
class HierarchicalClusterNode;
class HierarchicalCluster
{
	public:
	HierarchicalCluster();
	~HierarchicalCluster();
	int setOutputDir(const char*);
	int cluster(map<int,map<string,int>*>& modules,map<string,HierarchicalClusterNode*>& attribs,double t,map<int,HierarchicalClusterNode*>& attribs_integer,map<string,int>&);
	HierarchicalClusterNode* getRoot();
	double** getDist();
	private:
	int estimatePairwiseDist(map<int,HierarchicalClusterNode*>&,map<string,int>&);
	int mergePairs(map<int,HierarchicalClusterNode*>&);
	int mergePairs_LazyDelete(map<int,HierarchicalClusterNode*>&);
	int clearNeighborsFromList(HierarchicalClusterNode* node, map<string,int>& deleteMe);
	int generateModules(map<int,HierarchicalClusterNode*>& currNodeSet,map<int,map<string,int>*>& modules,map<string,HierarchicalClusterNode*>&);
	int calculatePercentVarianceExplained(map<int,map<string,int>*>& modules,map<string,HierarchicalClusterNode*>& origAttrib);
	int calculateSilhouetteIndex(map<int,map<string,int>*>& modules,map<string,HierarchicalClusterNode*>& origAttrib);
	int populateMembers(map<string,int>* members,HierarchicalClusterNode* node);
	int findRoot(HierarchicalClusterNode*n);
	double threshold;
	map<string,HierarchicalClusterNode*> backup;
	char outputDir[1024];
	HierarchicalClusterNode* root;
	int globalNodeID;
	int globalPairID;
	double** distvalues;
	int** neighbors;
	int* visited;
	int treenodecnt;
	struct Pair
	{
		int node1;
		int node2;
		double value;
	};
	class ComparePairMin
	{
		public:
		bool operator() (Pair* p1, Pair* p2)
		{
			if(p1->value > p2->value)
			{
				return true;
			}
			return false;
		}
	};
	class ComparePairMax
	{
		public:
		bool operator() (Pair* p1, Pair* p2)
		{
			if(p1->value < p2->value)
			{
				return true;
			}
			return false;
		}
	};
	priority_queue<HierarchicalCluster::Pair*,vector<HierarchicalCluster::Pair*>, HierarchicalCluster::ComparePairMax> myqueue;
	
};
#endif
