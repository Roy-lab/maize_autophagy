#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include "HierarchicalClusterNode.H"
#include "HierarchicalCluster.H"
#include "OptimalLeafOrder.H"
#include "Framework.H"

Framework::Framework()
{
}

Framework::~Framework()
{
}

int 
Framework::readDataMatrix(const char* aFName)
{
	ifstream inFile(aFName);
	char* buffer=NULL;
	int bufflen=0;
	string strbuff;
	int gid=0;
	while(inFile.good())
	{
		getline(inFile,strbuff);
		if(strbuff.length()<=0)
		{
			continue;
		}
		if(bufflen<=strbuff.length())
		{
			if(buffer!=NULL)
			{
				delete[] buffer;
			}
			bufflen=strbuff.length()+1;
			buffer=new char[bufflen];
		}
		strcpy(buffer,strbuff.c_str());
		if((strstr(buffer,"Gene")!=NULL) || strstr(buffer,"ORF")!=NULL)
		{
			readColumns(buffer);
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;	
		string geneName;
		HierarchicalClusterNode* node=new HierarchicalClusterNode;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				geneName.append(tok);
			}
			else
			{
				double v=atof(tok);
				if(v>0)
				{
					node->attrib[tokCnt-1]=1;
				}
				else
				{
					node->attrib[tokCnt-1]=0;
				}
				node->attrib[tokCnt-1]=v;
				node->expr.push_back(v);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		node->nodeName.append(geneName.c_str());
		nodeSet[geneName]=node;
		backup[geneName]=node;
		nameIDMap[geneName]=gid;
		node->size=1;
		gid++;
	}
	cout <<"Read " << nodeSet.size() << " points" << endl;
	inFile.close();
	return 0;
}

int
Framework::readDataPair(const char* aFName)
{
	ifstream inFile(aFName);
	char* buffer=NULL;
	int bufflen=0;
	string strbuff;
	int gid=1;
	while(inFile.good())
	{
		getline(inFile,strbuff);
		if(strbuff.length()<=0)
		{
			continue;
		}
		if(bufflen<=strbuff.length())
		{
			if(buffer!=NULL)
			{
				delete[] buffer;
			}
			bufflen=strbuff.length()+1;
			buffer=new char[bufflen];
		}
		strcpy(buffer,strbuff.c_str());
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;	
		string geneName;
		string attribName;
		double val=0;
		HierarchicalClusterNode* node=NULL;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				geneName.append(tok);
			}
			else if(tokCnt==1)
			{
				attribName.append(tok);
			}
			else if(tokCnt==2)
			{
				val=atof(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		int attribID=-1;
		if(nameIDMap.find(attribName)==nameIDMap.end())
		{
			attribID=nameIDMap.size();
			nameIDMap[attribName]=attribID;
		}
		else
		{
			attribID=nameIDMap[attribName];
		}
		if(nodeSet.find(geneName)==nodeSet.end())
		{
			node=new HierarchicalClusterNode;
			node->nodeName.append(geneName.c_str());
			nodeSet[geneName]=node;
			backup[geneName]=node;
			nameIDMap[geneName]=gid;
			gid++;
			node->size=1;
		}
		else
		{
			node=nodeSet[geneName];
		}
		node->attrib[attribID]=val;
	}
	cout <<"Read " << nodeSet.size() << " points" << endl;
	for(map<string,HierarchicalClusterNode*>::iterator hIter=nodeSet.begin();hIter!=nodeSet.end();hIter++)
	{
		cout << nameIDMap[hIter->first]<<"\t"<< hIter->first;
		HierarchicalClusterNode* hc=hIter->second;
		for(map<int,double>::iterator aIter=hc->attrib.begin();aIter!=hc->attrib.end();aIter++)
		{
			cout <<"\t" << aIter->second;
		}
		cout << endl;
	
	}
	inFile.close();
	return 0;
}




int 
Framework::reorder(const char* aFSuff,double threshold)
{
	map<int,map<string,int>*> modules;
	map<int,HierarchicalClusterNode*> attribs;
	cluster.cluster(modules,nodeSet,threshold,attribs,attribNameIDMap);
	double** dist=cluster.getDist();
	char aFName[1024];
	sprintf(aFName,"%s_assign.txt",aFSuff);
	ofstream oFile(aFName);
	sprintf(aFName,"%s_geneset.txt",aFSuff);
	ofstream gFile(aFName);
	int c=0;
	int atcnt=0;
	olo.setDist(dist);
	for(map<int,HierarchicalClusterNode*>::iterator aIter=attribs.begin();aIter!=attribs.end();aIter++)
	{
		olo.setHierarchicalClusterNode(aIter->second);
		vector<string> ordering;
		olo.reorder(ordering);
		for(int i=0;i<ordering.size();i++)
		{
			HierarchicalClusterNode* hc=backup[ordering[i]];
			oFile<<ordering[i];
			gFile <<ordering[i]<< "\t"<< c << endl;
			
			for(map<int,double>::iterator aIter=hc->attrib.begin();aIter!=hc->attrib.end();aIter++)
			{
				oFile<<"\t" << aIter->second;
			}
			oFile  <<endl;
			atcnt=hc->attrib.size();	
		}		
		//gFile << endl;
		oFile <<"Dummy"<< c;
		for(int i=0;i<atcnt;i++)
		{
			oFile <<"\t-999";
		}
		oFile << endl;

		
		c++;
	}
	oFile.close();
	gFile.close();
	return 0;
}


int
Framework::readColumns(char* buffer)
{
	int tokCnt=0;
	char* tok=strtok(buffer,"\t");
	while(tok!=NULL)
	{
		if(tokCnt>0)
		{
			string colName(tok);
			attribNameIDMap[colName]=tokCnt-1;
			attribIDNameMap[tokCnt-1]=colName;
		}
		tok=strtok(NULL,"\t");
		tokCnt++;
	}
	return 0;
}

int
main(int argc, const char** argv)
{
	if(argc!=5)
	{
		cout <<"Usage: getOrdering data [list|pair] outfile threshold" << endl;
		return 0;
	}
	Framework fw;
	if(strcmp(argv[2],"matrix")==0)
	{
		fw.readDataMatrix(argv[1]);
	}
	else
	{
		fw.readDataPair(argv[1]);
	}
	fw.reorder(argv[3],atof(argv[4]));
	return 0;
}

