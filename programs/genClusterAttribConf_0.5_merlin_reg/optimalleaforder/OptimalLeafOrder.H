#ifndef _OPTIMALLEAFORDER_
#define _OPTIMALLEAFORDER_
#include <map>
#include <string>
#include <vector>
using namespace std;
class HierarchicalClusterNode;
class OptimalLeafOrder
{
	public:
	OptimalLeafOrder();
	~OptimalLeafOrder();
	int setHierarchicalClusterNode(HierarchicalClusterNode*);
	int reorder(vector<string>&);
	private:
	typedef struct
	{
		string leftextreme;
		string rightextreme;
	}Pair;
	struct Node;
	struct Node
	{
		string nodeName;
		map<string,map<string,double>*> childOrder;
		map<string,Pair*> extremes;
		Node* left;
		Node* right;
		Node* parent;
	};
	Node* root;
	map<string,Node*> minDist;
	int getordering(Node*,string&,string&,vector<string>&);
	int populateDist(HierarchicalClusterNode* node);
	double getSim(string&,string&);
	int makeTree(HierarchicalClusterNode*,Node*);
	int reorder(Node* n);
	map<string,double> distance;
	HierarchicalClusterNode* hcRoot;
};
#endif
