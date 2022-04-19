#ifndef _HIERARCHICAL_CLUSTER
#define _HIERARCHICAL_CLUSTER
#include <map>
#include <string>
#include <list>
#include <vector>
using namespace std;
class HierarchicalClusterNode;
class Heap;
class HierarchicalCluster
{
	public:
	HierarchicalCluster();
	~HierarchicalCluster();
	int setOutputDir(const char*);
	int cluster(map<int,map<string,int>*>& modules,map<string,HierarchicalClusterNode*>& attribs,double);
	int clusterExp(map<int,map<string,int>*>& modules,map<string,HierarchicalClusterNode*>& attribs,double);
	HierarchicalClusterNode* getRoot();
	private:
	int estimatePairwiseDist(map<string,HierarchicalClusterNode*>&);
	int estimatePairwiseDistExp(map<string,HierarchicalClusterNode*>&);
	int mergePairs(map<string,HierarchicalClusterNode*>&);
	int mergePairs_Efficient(map<string,HierarchicalClusterNode*>&);
	int clearNeighborsFromList(HierarchicalClusterNode* node, map<string,int>& deleteMe);
	int generateModules(map<string,HierarchicalClusterNode*>& currNodeSet,map<int,map<string,int>*>& modules,map<string,HierarchicalClusterNode*>&);
	int calculatePercentVarianceExplained(map<int,map<string,int>*>& modules,map<string,HierarchicalClusterNode*>& origAttrib);
	int calculateSilhouetteIndex(map<int,map<string,int>*>& modules,map<string,HierarchicalClusterNode*>& origAttrib);
	int populateMembers(map<string,int>* members,HierarchicalClusterNode* node);
	int findRoot(HierarchicalClusterNode*n);
	double threshold;
	Heap heap;
	map<string,Heap*> pairPtrMap;
	map<string,HierarchicalClusterNode*> backup;
	char outputDir[1024];
	HierarchicalClusterNode* root;
};
#endif
