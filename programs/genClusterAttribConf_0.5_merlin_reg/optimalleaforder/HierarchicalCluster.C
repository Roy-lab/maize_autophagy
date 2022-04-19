#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include <time.h>
#include "Distance.H"
#include "HierarchicalClusterNode.H"
#include "Heap.H"
#include "HierarchicalCluster.H"
int sortfunc(const void* first, const void* second);
int* sortingind=NULL;
double* sorteddist=NULL;


HierarchicalCluster::HierarchicalCluster()
{
	root=NULL;
}

HierarchicalCluster::~HierarchicalCluster()
{
}


int
HierarchicalCluster::setOutputDir(const char* aFName)
{
	strcpy(outputDir,aFName);
	return 0;
}

int 
HierarchicalCluster::clusterExp(map<int,map<string,int>*>& modules,map<string,HierarchicalClusterNode*>& attribs,double t)
{
	threshold=t;
	estimatePairwiseDistExp(attribs);
	for(map<string,HierarchicalClusterNode*>::iterator aIter=attribs.begin();aIter!=attribs.end();aIter++)
	{
		backup[aIter->first]=aIter->second;
	}
	
	bool keepMerging=true;
	
	while(keepMerging && attribs.size()>1)
	{
		struct timeval begintime;
		struct timeval endtime;
		gettimeofday(&begintime,NULL);
		//keepMerging=mergePairs(attribs);
		keepMerging=mergePairs_Efficient(attribs);
		gettimeofday(&endtime,NULL);
		//cout << "Time elapsed " << endtime.tv_sec-begintime.tv_sec<< " seconds and " << endtime.tv_usec-begintime.tv_usec << " micro secs" << endl;
	}
	generateModules(attribs,modules,backup);
	calculatePercentVarianceExplained(modules,backup);
	calculateSilhouetteIndex(modules,backup);
	/*for(map<string,HierarchicalClusterNode*>::iterator aIter=backup.begin();aIter!=backup.end();aIter++)
	{
		map<string,HierarchicalClusterNode*>::iterator cIter=attribs.find(aIter->first);
		if(cIter!=attribs.end())
		{
			attribs.erase(cIter);
		}
		delete aIter->second;
	}
	for(map<string,HierarchicalClusterNode*>::iterator aIter=attribs.begin();aIter!=attribs.end();aIter++)
	{
		delete aIter->second;
	}
	attribs.clear();
	backup.clear();
	pairPtrMap.clear();*/
	return 0;
}

int 
HierarchicalCluster::cluster(map<int,map<string,int>*>& modules,map<string,HierarchicalClusterNode*>& attribs,double t)
{
	threshold=t;
	estimatePairwiseDist(attribs);
	for(map<string,HierarchicalClusterNode*>::iterator aIter=attribs.begin();aIter!=attribs.end();aIter++)
	{
		backup[aIter->first]=aIter->second;
	}
	
	bool keepMerging=true;
	
	while(keepMerging && attribs.size()>1)
	{
		struct timeval begintime;
		struct timeval endtime;
		gettimeofday(&begintime,NULL);
		//keepMerging=mergePairs(attribs);
		keepMerging=mergePairs_Efficient(attribs);
		gettimeofday(&endtime,NULL);
		//cout << "Time elapsed " << endtime.tv_sec-begintime.tv_sec<< " seconds and " << endtime.tv_usec-begintime.tv_usec << " micro secs" << endl;
	}
	generateModules(attribs,modules,backup);
	calculatePercentVarianceExplained(modules,backup);
	calculateSilhouetteIndex(modules,backup);
	/*for(map<string,HierarchicalClusterNode*>::iterator aIter=backup.begin();aIter!=backup.end();aIter++)
	{
		map<string,HierarchicalClusterNode*>::iterator cIter=attribs.find(aIter->first);
		if(cIter!=attribs.end())
		{
			attribs.erase(cIter);
		}
		delete aIter->second;
	}
	for(map<string,HierarchicalClusterNode*>::iterator aIter=attribs.begin();aIter!=attribs.end();aIter++)
	{
		delete aIter->second;
	}
	attribs.clear();
	backup.clear();
	pairPtrMap.clear();*/
	return 0;
}

int
HierarchicalCluster::estimatePairwiseDistExp(map<string,HierarchicalClusterNode*>& currNodeSet)
{
	Distance d;
	map<int,string> nodeIDName;
	map<int,double> nodeIDDist;
	map<string,double> nodeNameDist;
	for(map<string,HierarchicalClusterNode*>::iterator nIter=currNodeSet.begin();nIter!=currNodeSet.end();nIter++)
	{
		HierarchicalClusterNode* hcNode1=nIter->second;
		map<string,HierarchicalClusterNode*>::iterator mIter=nIter;
		mIter++;
		for(;mIter!=currNodeSet.end();mIter++)
		{
			HierarchicalClusterNode* hcNode2=mIter->second;
			double rdist=0;
			double ccdist=0;
			double sharedSign=0;
			double den1=0;
			double den2=0;
			double dist=d.computeCC(hcNode1->expr,hcNode2->expr);
			dist=0.5*(1-dist);
			hcNode1->distToNeighbors[hcNode2->nodeName]=dist;
			hcNode2->distToNeighbors[hcNode1->nodeName]=dist;
			string nodePairName;
			char buffer[1024];
			sprintf(buffer,"%s-%s",hcNode1->nodeName.c_str(),hcNode2->nodeName.c_str());
			Heap* ptr=heap.insertToHeapNoHeapify(hcNode1->nodeName,hcNode2->nodeName,dist);
			nodePairName.append(buffer);
			pairPtrMap[nodePairName]=ptr;
			int size=nodeIDName.size();
			nodeIDName[size]=nodePairName;
			//nodeIDDist[size]=dist;
			nodeNameDist[nodePairName]=dist;
		}
	}
	sorteddist=new double[nodeIDName.size()];
	sortingind=new int[nodeIDName.size()];
	for(map<int,string>::iterator nIter=nodeIDName.begin();nIter!=nodeIDName.end();nIter++)
	{
		sortingind[nIter->first]=nIter->first;
		sorteddist[nIter->first]=nodeNameDist[nIter->second];
	}
	if(!heap.checkHeap())
	{
		cout <<"Heap violated!" << endl;
	}
	if(!heap.checkPointers(heap.getRoot()))
	{
		cout <<"Heap pointers off" << endl;	
	}
	//heap.showHeap();
	/*for(int i=0;i<nodeIDName.size();i++)
	{
		int sid=sortingind[i];
		string& pairName=nodeIDName[sid];
		Pair* p=new Pair;
		int pos=pairName.find("-");
		p->node1.append(pairName.substr(0,pos));
		p->node2.append(pairName.substr(pos+1,pairName.length()-pos));
		p->dist=nodeNameDist[pairName];
		sortedListNodePair.push_back(p);
		sortedVectorNodePair.push_back(p);
		n1->neighborIter[pairName]=sortedVectorNodePair.size()-1;
		HierarchicalClusterNode* n2=currNodeSet[p->node2];
		n2->neighborIter[pairName]=sortedVectorNodePair.size()-1;
	}*/
	delete[] sortingind;
	delete[] sorteddist;
	nodeIDName.clear();
	nodeNameDist.clear();
	/*for(list<Pair*>::iterator lIter=sortedListNodePair.begin();lIter!=sortedListNodePair.end();lIter++)
	{
		cout << (*lIter)->node1 << "-" << (*lIter)->node2 <<"\t" << (*lIter)->dist << endl;
	}*/
	//make_heap(sortedVectorNodePair.begin(),sortedVectorNodePair.end());
	return 0;
}

int
HierarchicalCluster::estimatePairwiseDist(map<string,HierarchicalClusterNode*>& currNodeSet)
{
	Distance d;
	map<int,string> nodeIDName;
	map<int,double> nodeIDDist;
	map<string,double> nodeNameDist;
	for(map<string,HierarchicalClusterNode*>::iterator nIter=currNodeSet.begin();nIter!=currNodeSet.end();nIter++)
	{
		HierarchicalClusterNode* hcNode1=nIter->second;
		map<string,HierarchicalClusterNode*>::iterator mIter=nIter;
		mIter++;
		for(;mIter!=currNodeSet.end();mIter++)
		{
			HierarchicalClusterNode* hcNode2=mIter->second;
			double rdist=0;
			double ccdist=0;
			double sharedSign=0;
			double den1=0;
			double den2=0;
			for(map<int,double>::iterator aIter=hcNode1->attrib.begin();aIter!=hcNode1->attrib.end();aIter++)
			{
				//rdist=rdist+(fabs(aIter->second-hcNode2->attrib[aIter->first]));
				double othernodeval=0;
				if(hcNode2->attrib.find(aIter->first)==hcNode2->attrib.end())
				{	
					rdist=rdist+1;
				}
			}
			for(map<int,double>::iterator aIter=hcNode2->attrib.begin();aIter!=hcNode2->attrib.end();aIter++)
			{
				if(hcNode1->attrib.find(aIter->first)!=hcNode1->attrib.end())
				{
					continue;
				}
				rdist=rdist+1;
			}
			//double dist=d.computeCC(hcNode1->expr,hcNode2->expr);
			//dist=0.5*(1-dist);
			double dist=rdist;
			hcNode1->distToNeighbors[hcNode2->nodeName]=dist;
			hcNode2->distToNeighbors[hcNode1->nodeName]=dist;
			string nodePairName;
			char buffer[1024];
			sprintf(buffer,"%s-%s",hcNode1->nodeName.c_str(),hcNode2->nodeName.c_str());
			Heap* ptr=heap.insertToHeapNoHeapify(hcNode1->nodeName,hcNode2->nodeName,dist);
			nodePairName.append(buffer);
			pairPtrMap[nodePairName]=ptr;
			int size=nodeIDName.size();
			nodeIDName[size]=nodePairName;
			//nodeIDDist[size]=dist;
			nodeNameDist[nodePairName]=dist;
		}
	}
	sorteddist=new double[nodeIDName.size()];
	sortingind=new int[nodeIDName.size()];
	for(map<int,string>::iterator nIter=nodeIDName.begin();nIter!=nodeIDName.end();nIter++)
	{
		sortingind[nIter->first]=nIter->first;
		sorteddist[nIter->first]=nodeNameDist[nIter->second];
	}
	if(!heap.checkHeap())
	{
		cout <<"Heap violated!" << endl;
	}
	if(!heap.checkPointers(heap.getRoot()))
	{
		cout <<"Heap pointers off" << endl;	
	}
	//heap.showHeap();
	/*for(int i=0;i<nodeIDName.size();i++)
	{
		int sid=sortingind[i];
		string& pairName=nodeIDName[sid];
		Pair* p=new Pair;
		int pos=pairName.find("-");
		p->node1.append(pairName.substr(0,pos));
		p->node2.append(pairName.substr(pos+1,pairName.length()-pos));
		p->dist=nodeNameDist[pairName];
		sortedListNodePair.push_back(p);
		sortedVectorNodePair.push_back(p);
		n1->neighborIter[pairName]=sortedVectorNodePair.size()-1;
		HierarchicalClusterNode* n2=currNodeSet[p->node2];
		n2->neighborIter[pairName]=sortedVectorNodePair.size()-1;
	}*/
	delete[] sortingind;
	delete[] sorteddist;
	nodeIDName.clear();
	nodeNameDist.clear();
	/*for(list<Pair*>::iterator lIter=sortedListNodePair.begin();lIter!=sortedListNodePair.end();lIter++)
	{
		cout << (*lIter)->node1 << "-" << (*lIter)->node2 <<"\t" << (*lIter)->dist << endl;
	}*/
	//make_heap(sortedVectorNodePair.begin(),sortedVectorNodePair.end());
	return 0;
}

int
HierarchicalCluster::mergePairs(map<string,HierarchicalClusterNode*>& currNodeSet)
{
	double maxSum=0;
	string node1;
	string node2;
	double theoreticalMax=1;
	double minDist=100000;
	for(map<string,HierarchicalClusterNode*>::iterator nIter=currNodeSet.begin();nIter!=currNodeSet.end();nIter++)
	{
		HierarchicalClusterNode* hcNode1=nIter->second;
		if(hcNode1->status==1)
		{
			continue;
		}
		map<string,HierarchicalClusterNode*>::iterator mIter=nIter;
		mIter++;
		for(;mIter!=currNodeSet.end();mIter++)
		{
			HierarchicalClusterNode* hcNode2=mIter->second;
			if(hcNode2->status==1)
			{
				continue;
			}
			/*double sim=0;
			for(map<int,double>::iterator aIter=hcNode1->attrib.begin();aIter!=hcNode1->attrib.end();aIter++)
			{
				if(hcNode2->attrib.find(aIter->first)!=hcNode2->attrib.end())
				{
					sim=sim+1;
				}
			}
			sim=sim/(double(hcNode1->attrib.size()+hcNode2->attrib.size()-sim));*/
			if(hcNode1->distToNeighbors.find(mIter->first)==hcNode1->distToNeighbors.end())
			{
				cout <<"No distance of " << mIter->first << " in " << hcNode1->nodeName << endl; 
				exit(0);
			}
			double dist=hcNode1->distToNeighbors[mIter->first];
			if(dist<minDist)
			{
				if(node1.length()>0)
				{	
					node1.clear();
				}
				if(node2.length()>0)
				{
					node2.clear();
				}
				node1.append(nIter->first.c_str());
				node2.append(mIter->first.c_str());
				minDist=dist;
				//if(maxSum>=theoreticalMax)
				//{
				//	break;
				//}
			}
		}
		
		if(maxSum>=theoreticalMax)
		{
			break;
		}
	}
	if(minDist>=threshold)
	{
	//	return false;
	}
	if(currNodeSet.find(node1)==currNodeSet.end())
	{
		cout <<"node1 " << node1 << " not found " << endl;
		exit(0);
	}
	if(currNodeSet.find(node2)==currNodeSet.end())
	{
		cout <<"node2 " << node2 << " not found " << endl;
		exit(0);
	}

	HierarchicalClusterNode* c1=currNodeSet[node1];
	HierarchicalClusterNode* c2=currNodeSet[node2];
	//c1->status=1;
	//c2->status=1;
	HierarchicalClusterNode* c12=new HierarchicalClusterNode;
	/*for(map<int,double>::iterator aIter=c1->attrib.begin();aIter!=c1->attrib.end();aIter++)
	{
		if(c2->attrib.find(aIter->first)==c2->attrib.end())
		{
			c12->attrib[aIter->first]=aIter->second;
			continue;
		}
		if((aIter->second*c1->attrib[aIter->first]) < 0)
		{
			c12->attrib[aIter->first]=(aIter->second+c1->attrib[aIter->first])/2;
			continue;
		}	
	//	c12->attrib[aIter->first]=aIter->second;
		c12->attrib[aIter->first]=(aIter->second+c1->attrib[aIter->first])/2;
	}
	for(map<int,double>::iterator aIter=c2->attrib.begin();aIter!=c2->attrib.end();aIter++)
	{
		if(c1->attrib.find(aIter->first)==c1->attrib.end())
		{
			c12->attrib[aIter->first]=aIter->second;
		}
	}
	if(c12->attrib.size()<threshold)
	{
		cout << "Merge not possible for " << node1 << " " << node2 << " " << minDist << endl;
		//don't merge move to the next candidate.
		return true;
		delete c12;
	}
	for(int i=0;i<c1->expr.size();i++)
	{
		double eval=(c1->expr[i]+c2->expr[i])/2;
		c12->expr.push_back(eval);
	}*/
	cout << "Merging " << node1 << " " << node2  << " dist=" << minDist << endl;
	c12->left=c1;
	c12->right=c2;
	c1->parent=c12;
	c2->parent=c12;
	c12->nodeName.append(node1);
	c12->nodeName.append("-");
	c12->nodeName.append(node2);
	/*for(map<int,int>::iterator aIter=c2->attrib.begin();aIter!=c2->attrib.end();aIter++)
	{
		c12->attrib[aIter->first]=aIter->second;
	}*/
	map<string,HierarchicalClusterNode*>::iterator hIter1=currNodeSet.find(node1);
	if(hIter1==currNodeSet.end())
	{
		cout << "Something funny happened. Can't find node to delete!" << endl;
		exit(0);
	}
	//delete hIter1->second;
	currNodeSet.erase(hIter1);
	map<string,HierarchicalClusterNode*>::iterator hIter2=currNodeSet.find(node2);
	if(hIter2==currNodeSet.end())
	{
		cout << "Something funny happened. Can't find node to delete!" << endl;
		exit(0);
	}
	//delete hIter2->second;
	currNodeSet.erase(hIter2);
	Distance cc;
	for(map<string,HierarchicalClusterNode*>::iterator nIter=currNodeSet.begin();nIter!=currNodeSet.end();nIter++)
	{
		HierarchicalClusterNode* n1=nIter->second;
		if(n1->status==1)
		{
			continue;
		}
		double d1=c1->distToNeighbors[n1->nodeName];
		double d2=c2->distToNeighbors[n1->nodeName];
		double dkm_rdist=((c1->size*d1) + (c2->size*d2))/((double)(c1->size+ c2->size));
		double sim=0;
		double den1=0;
	/*	for(map<int,double>::iterator aIter=c12->attrib.begin();aIter!=c12->attrib.end();aIter++)
		{
			den1=den1+fabs(aIter->second);
			if(n1->attrib.find(aIter->first)==n1->attrib.end())
			{
				continue;
			}
			if((aIter->second*n1->attrib[aIter->first])>=0)
			{
				sim=sim+((fabs(aIter->second)+fabs(n1->attrib[aIter->first]))/2.0);
				//sim=sim+1;
			}
		}
		double den2=0;
		for(map<int,double>::iterator aIter=n1->attrib.begin();aIter!=n1->attrib.end();aIter++)
		{
			den2=den2+fabs(aIter->second);
		}*/
		//dkm=1-(sim/((double) (c12->attrib.size()+n1->attrib.size()-sim)));
		//double dkm_rdist=1-(sim/((double) (den1+den2-sim)));
		//double dkm_ccdist=(1-cc.computeCC(c12->expr,n1->expr))/2.0;
		//dkm_rdist=dkm_ccdist;
		//double dist=0.5*(dkm_rdist+dkm_ccdist);
		double dist=dkm_rdist;
		c12->distToNeighbors[n1->nodeName]=dist;
		n1->distToNeighbors[c12->nodeName]=dist;
		//if(strcmp(n1->nodeName.c_str(),"YFR055W")==0)
		if(dist<threshold)
		{
			//cout <<n1->nodeName<<"\t" << c12->nodeName << " dist " << dist << " CC=" << dkm_ccdist << " Reg=" << dkm_rdist <<endl;
			//cout << "newdists " << dist << " CC=" << dkm_ccdist << " Reg=" << dkm_rdist <<endl;
		}
	}
	currNodeSet[c12->nodeName]=c12;
	c12->size=c1->size+c2->size;
	return true;
}


int
HierarchicalCluster::mergePairs_Efficient(map<string,HierarchicalClusterNode*>& currNodeSet)
{
	double maxSum=0;
	double minDist=100000;
	if(heap.empty())
	{
		return false;
	}
	Heap::Pair* p=heap.getMin();
	minDist=p->dist;
	if(p->dist>=threshold)
	{
//		return false;
	}
	cout << "Merging " << p->node1 << " " << p->node2  << " dist=" << minDist << endl;
	if((strcmp(p->node1.c_str(),"ADH2-BAR1-SPO11-FUS1-SIN3-RME1-SWI1")==0) || (strcmp(p->node2.c_str(),"ALPHA1-HXT7-MFALPHA2-SNF2_SWI1-ARG5-HXT6-MCM1-STA2-STA1-SAG1-STE6-PHO11-PHO5")==0))
	{
		cout << "Found CDC9" << endl;
		Heap::Pair* p1=heap.getMin();
	}
	
	if(currNodeSet.find(p->node1)==currNodeSet.end())
	{
		cout <<"node1 " << p->node1 << " not found in pair " << p->node1<<"-" <<p->node2 << endl;
		exit(0);
	}
	if(currNodeSet.find(p->node2)==currNodeSet.end())
	{
		cout <<"node2 " << p->node2 << " not found in pair " << p->node1 <<"-" << p->node2<<  endl;
		exit(0);
	}
	HierarchicalClusterNode* c1=currNodeSet[p->node1];
	HierarchicalClusterNode* c2=currNodeSet[p->node2];
	HierarchicalClusterNode* c12=new HierarchicalClusterNode;
	c12->left=c1;
	c12->right=c2;
	c1->parent=c12;
	c2->parent=c12;
	c12->nodeName.append(p->node1);
	c12->nodeName.append("-");
	c12->nodeName.append(p->node2);
	map<string,HierarchicalClusterNode*>::iterator hIter1=currNodeSet.find(c1->nodeName);
	if(hIter1==currNodeSet.end())
	{
		cout << "Something funny happened. Can't find node to delete!" << endl;
		exit(0);
	}
	HierarchicalClusterNode* n1=hIter1->second;
	for(map<string,double>::iterator sIter=n1->distToNeighbors.begin();sIter!=n1->distToNeighbors.end();sIter++)
	{
		string key;
		if(strcmp(sIter->first.c_str(),n1->nodeName.c_str())<0)
		{
			key.append(sIter->first);
			key.append("-");
			key.append(n1->nodeName);
		}	
		else
		{
			key.append(n1->nodeName);
			key.append("-");
			key.append(sIter->first);
		}
		//cout <<"Deleting " << key << endl;
		map<string,Heap*>::iterator heapIter=pairPtrMap.find(key);
		if(heapIter==pairPtrMap.end())
		{
			continue;
		}
		Heap* todel=heapIter->second;
		//heap.deleteFromHeap(todel);
		heap.deleteFromHeap_getLeaf(todel);
		/*if(!heap.checkPointers())
		{
			cout <<"Heap pointers off" << endl;	
		}*/
		pairPtrMap.erase(heapIter);
	}
	//clearNeighborsFromList(hIter1->second,deleteMe);
	//heap.deleteFromHeap(heap.getRoot());
	map<string,HierarchicalClusterNode*>::iterator hIter2=currNodeSet.find(c2->nodeName);
	if(hIter2==currNodeSet.end())
	{
		cout << "Something funny happened. Can't find node to delete!" << endl;
		exit(0);
	}
	//clearNeighborsFromList(hIter2->second,deleteMe);
	HierarchicalClusterNode* n2=hIter2->second;
	for(map<string,double>::iterator sIter=n2->distToNeighbors.begin();sIter!=n2->distToNeighbors.end();sIter++)
	{
		string key;
		if(strcmp(sIter->first.c_str(),n2->nodeName.c_str())<0)
		{
			key.append(sIter->first);
			key.append("-");
			key.append(n2->nodeName);
		}	
		else
		{
			key.append(n2->nodeName);
			key.append("-");
			key.append(sIter->first);
		}
		if(strcmp(n2->nodeName.c_str(),"CDC9")==0)
		{
			cout <<"Found CDC9 other neighbor " << key << endl;
			if(strcmp(key.c_str(),"CDC9-RAD27")==0)
			{
				cout << "Stop here" << endl;
			}
		}
		map<string,Heap*>::iterator heapIter=pairPtrMap.find(key);
		if(heapIter==pairPtrMap.end())
		{
			continue;
		}
		Heap* todel=heapIter->second;
		//heap.deleteFromHeap(todel);
		heap.deleteFromHeap_getLeaf(todel);
		/*if(!heap.checkPointers())
		{
			cout <<"Heap pointers off" << endl;	
		}*/
		pairPtrMap.erase(heapIter);
	}
	/*list<Pair*>::iterator dIter=sortedListNodePair.begin();
	cout <<"To delete " << deleteMe.size()  << " from " << sortedListNodePair.size()<< endl;
	while(dIter!=sortedListNodePair.end())
	{	
		Pair* p=(*dIter);
		string key(p->node1);
		key.append("-");
		key.append(p->node2);
		list<Pair*>::iterator nIter=dIter;
		nIter++;
		if(deleteMe.find(key)!=deleteMe.end())
		{
			sortedListNodePair.erase(dIter);
		}
		dIter=nIter;
	}
	deleteMe.clear();
	cout <<"Post deleting " << sortedListNodePair.size()<< " nodes" << endl;*/
	if(strcmp(hIter1->first.c_str(),"CDC9")==0)
	{
		cout <<"Deleting CDC9" << endl;
	}
	if(strcmp(hIter2->first.c_str(),"CDC9")==0)
	{
		cout <<"Deleting CDC9" << endl;
	}
	currNodeSet.erase(hIter1);
	currNodeSet.erase(hIter2);
	//cout <<"Prior to adding new pairs" << endl;
	//heap.showHeap();
	for(map<string,HierarchicalClusterNode*>::iterator nIter=currNodeSet.begin();nIter!=currNodeSet.end();nIter++)
	{
		HierarchicalClusterNode* n1=nIter->second;
		if(n1->status==1)
		{
			continue;
		}
		double d1=c1->distToNeighbors[n1->nodeName];
		double d2=c2->distToNeighbors[n1->nodeName];
		double dkm_rdist=((c1->size*d1) + (c2->size*d2))/((double)(c1->size+ c2->size));
		double dist=dkm_rdist;
		c12->distToNeighbors[n1->nodeName]=dist;
		n1->distToNeighbors[c12->nodeName]=dist;
		if(dist>threshold)
		{
			//continue;
		}
		//Pair* newp=new Pair;
		string key;
		if(strcmp(c12->nodeName.c_str(),n1->nodeName.c_str())<0)
		{
			Heap* newp=heap.insertToHeapNoHeapify(c12->nodeName,n1->nodeName,dist);
			key.append(c12->nodeName);
			key.append("-");
			key.append(n1->nodeName);
			pairPtrMap[key]=newp;
		}
		else
		{
			Heap* newp=heap.insertToHeapNoHeapify(n1->nodeName,c12->nodeName,dist);
			key.append(n1->nodeName);
			key.append("-");
			key.append(c12->nodeName);
			pairPtrMap[key]=newp;
		}
		if(strcmp(key.c_str(),"CDC9-RAD27")==0)
		{
			cout << "Adding CDC9-RAD27 again!" << endl;
		}
		/*list<Pair*>::iterator lIter=sortedListNodePair.begin();
		bool greater=false;
		while(lIter!=sortedListNodePair.end() && !greater)
		{
			if((*lIter)->dist>=dist)
			{
				greater=true;
			}
			else
			{
				lIter++;
			}
		}
		if(greater)
		{
			sortedListNodePair.insert(lIter,newp);
		}
		else
		{
			sortedListNodePair.push_back(newp);
		}*/
	}
	/*cout <<"Post adding new nodes: " << sortedListNodePair.size() << endl;
	for(list<Pair*>::iterator lIter=sortedListNodePair.begin();lIter!=sortedListNodePair.end();lIter++)
	{
		cout << (*lIter)->node1 << "-" << (*lIter)->node2 <<"\t" << (*lIter)->dist << endl;
	}*/
	currNodeSet[c12->nodeName]=c12;
	c12->size=c1->size+c2->size;
	
	if(c1->left!=NULL || c1->right!=NULL)
	{
		backup[c1->nodeName]=c1;
	}
	if(c2->left!=NULL || c2->right!=NULL)
	{
		backup[c2->nodeName]=c2;
	}
	//cout <<"After adding pairs" << endl;
	//heap.showHeap();
	return true;
}

int
HierarchicalCluster::clearNeighborsFromList(HierarchicalClusterNode* node, map<string,int>& deleteMe)
{
	for(map<string,double>::iterator aIter=node->distToNeighbors.begin();aIter!=node->distToNeighbors.end();aIter++)
	{
		string key;
		if(strcmp(aIter->first.c_str(),node->nodeName.c_str())<0)
		{
			key.append(aIter->first);
			key.append("-");
			key.append(node->nodeName.c_str());
		}	
		else
		{
			key.append(node->nodeName.c_str());
			key.append("-");
			key.append(aIter->first);
		}
		deleteMe[key]=0;
	}
	node->attrib.clear();
	return 0;
}

int
HierarchicalCluster::generateModules(map<string,HierarchicalClusterNode*>& currNodeSet,map<int,map<string,int>*>& modules,map<string,HierarchicalClusterNode*>& origAttrib)
{
	int moduleCnt=modules.size();
	char outFName[1024];
	sprintf(outFName,"%s/clusters_pw.txt",outputDir);
	ofstream oFile(outFName);
	for(map<string,HierarchicalClusterNode*>::iterator cIter=currNodeSet.begin();cIter!=currNodeSet.end();cIter++)
	{
		HierarchicalClusterNode* node=cIter->second;
		if(node->status==1)
		{
			continue;
		}
		map<string,int>* moduleMembers=new map<string,int>;
		modules[moduleCnt]=moduleMembers;
		populateMembers(moduleMembers,node);
		for(map<string,int>::iterator aIter=moduleMembers->begin();aIter!=moduleMembers->end();aIter++)
		{
			HierarchicalClusterNode* childNode=origAttrib[aIter->first];
			oFile << aIter->first <<"||5\tModule||5\t" << moduleCnt <<"|1|"<< moduleCnt << endl;
			for(int i=0;i<childNode->expr.size();i++)
			{
				oFile <<aIter->first <<"||5\tExp"<< i << "||5\t" << childNode->expr[i] << "|2"<< endl; 
			}
			for(map<int,double>::iterator mIter=childNode->attrib.begin();mIter!=childNode->attrib.end();mIter++)
			{
				oFile << aIter->first <<"||5\t" << (mIter->first) << "||5\t1|2" << endl;
			}
		}
		oFile <<"|Spacer||"<<moduleCnt<<"|-" << endl;
		cout <<"Module: " << moduleCnt << "\tSize="<< moduleMembers->size() << endl;
	
		moduleCnt=moduleCnt+1;
	}
	oFile << endl;
	return 0;
}


int
HierarchicalCluster::calculatePercentVarianceExplained(map<int,map<string,int>*>& modules,map<string,HierarchicalClusterNode*>& origAttrib)
{
	double s_total=0;
	double s_err=0;
	map<int,double> globalMean;
	for(map<int,map<string,int>*>::iterator mIter=modules.begin();mIter!=modules.end();mIter++)
	{
		map<int,double> localMean;
		map<string,int>* geneset=mIter->second;
		for(map<string,int>::iterator gIter=geneset->begin();gIter!=geneset->end();gIter++)
		{
			HierarchicalClusterNode* n=origAttrib[gIter->first];
			for(int i=0;i<n->expr.size();i++)
			{
				if(localMean.find(i)==localMean.end())
				{
					localMean[i]=n->expr[i];
				}	
				else	
				{
					localMean[i]=localMean[i]+n->expr[i];
				}
			}
		}
		for(map<int,double>::iterator dIter=localMean.begin();dIter!=localMean.end();dIter++)
		{
			if(globalMean.find(dIter->first)==globalMean.end())
			{
				globalMean[dIter->first]=dIter->second;
			}
			else
			{
				globalMean[dIter->first]=globalMean[dIter->first]+dIter->second;
			}
			dIter->second=dIter->second/((double)geneset->size());
		}	
		double s_err_m=0;
		for(map<string,int>::iterator gIter=geneset->begin();gIter!=geneset->end();gIter++)
		{
			HierarchicalClusterNode* n=origAttrib[gIter->first];
			for(int i=0;i<n->expr.size();i++)
			{
				double diff=n->expr[i]-localMean[i];
				s_err_m=s_err_m+(diff*diff);
			}
		}
		s_err=s_err+s_err_m;
		localMean.clear();
	}
	for(map<int,double>::iterator dIter=globalMean.begin();dIter!=globalMean.end();dIter++)
	{
		dIter->second=dIter->second/((double)origAttrib.size());
	}
	for(map<int,map<string,int>*>::iterator mIter=modules.begin();mIter!=modules.end();mIter++)
	{
		map<string,int>* geneset=mIter->second;
		for(map<string,int>::iterator gIter=geneset->begin();gIter!=geneset->end();gIter++)
		{
			HierarchicalClusterNode* n=origAttrib[gIter->first];
			for(int i=0;i<n->expr.size();i++)
			{
				double diff=n->expr[i]-globalMean[i];
				s_total=s_total+(diff*diff);
			}
		}
	}
	double pcv=1.0-(s_err/s_total);
	cout <<"Percent variance explained " << pcv << endl;
	globalMean.clear();
	return 0;
}

int
HierarchicalCluster::calculateSilhouetteIndex(map<int,map<string,int>*>& modules,map<string,HierarchicalClusterNode*>& origAttrib)
{
	map<string,double> silhouette;
	int positive_s=0;
	double totals=0;
	for(map<int,map<string,int>*>::iterator mIter=modules.begin(); mIter!=modules.end();mIter++)
	{
		double module_s=0;
		map<string,int>* geneset=mIter->second;
		for(map<string,int>::iterator gIter=geneset->begin();gIter!=geneset->end();gIter++)
		{
			HierarchicalClusterNode* n=origAttrib[gIter->first];
			double a=0;
			for(map<string,int>::iterator hIter=geneset->begin();hIter!=geneset->end();hIter++)
			{
				if(gIter==hIter)
				{
					continue;
				}
				a=a+n->distToNeighbors_CC[hIter->first];
			}
			if(geneset->size()>1)
			{
				a=a/(geneset->size()-1);
			}
			double minb=100;
			for(map<int,map<string,int>*>::iterator nIter=modules.begin();nIter!=modules.end();nIter++)
			{
				double b=0;
				map<string,int>* geneset2=nIter->second;
				if(mIter==nIter)
				{
					continue;
				}
				for(map<string,int>::iterator hIter=geneset2->begin();hIter!=geneset2->end();hIter++)
				{
					b=b+n->distToNeighbors_CC[hIter->first];
				}
				b=b/geneset2->size();
				if(b<minb)
				{
					minb=b;
				}
			}
			double s=minb-a;
			if(s<0)
			{
				s=s/a;
			}
			else
			{
				s=s/minb;
				positive_s=positive_s+1;
			}
			silhouette[gIter->first]=s;
			module_s=module_s+s;
		}
		totals=totals+module_s;
		cout <<"Silhoutte for module " << mIter->first <<" " << module_s/geneset->size() << endl;
	}
	cout <<"Silhoutte Avg. " << totals/origAttrib.size() << " total positive " << positive_s << " of total " << origAttrib.size()<< endl;
	return 0;
}


int
HierarchicalCluster::populateMembers(map<string,int>* members,HierarchicalClusterNode* node)
{
	if(node->left==NULL && node->right==NULL)
	{
		(*members)[node->nodeName]=0;
	}
	else
	{
		if(node->left!=NULL)
		{
			populateMembers(members,node->left);
		}
		if(node->right!=NULL)
		{
			populateMembers(members,node->right);
		}
	}
	return 0;
}


HierarchicalClusterNode*
HierarchicalCluster::getRoot()
{
	if(root==NULL)
	{
		if(backup.size()==0)
		{
			cout <<"No nodes you fool!!" << endl;
			exit(0);
		}
		findRoot(backup.begin()->second);
	}
	return root;
}

int
HierarchicalCluster::findRoot(HierarchicalClusterNode* n)
{
	if(n->parent==NULL)
	{
		root=n;
	}
	else
	{
		findRoot(n->parent);
	}
	return 0;
}

int 
sortfunc(const void* first, const void* second)
{
	int ind1=*((int*)first);	
	int ind2=*((int*)second);
	double pval1=sorteddist[ind1];
	double pval2=sorteddist[ind2];
	int compstat=0;
	if(pval1<pval2)
	{
		compstat=-1;
	}
	else if(pval1>pval2)
	{
		compstat=1;
	}
	return compstat;
}
