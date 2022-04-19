#ifndef _FRAMEWORK_
#define _FRAMEWORK_

#include <map>
#include <string>
using namespace std;

class Framework
{
	public:
		Framework();
		~Framework();
		int readClusters(const char*);
		int showAdjacencyMap(const char*);
	private:
		int readCluster(const char*);
		map<string,map<string,double>*> adjacencymatrix;
		double clusterings;
};
#endif
