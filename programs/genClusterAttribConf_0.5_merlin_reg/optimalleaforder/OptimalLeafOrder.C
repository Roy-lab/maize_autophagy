#include <string.h>
#include <iostream>
#include <stdlib.h>
#include "HierarchicalClusterNode.H"
#include "OptimalLeafOrder.H"


OptimalLeafOrder::OptimalLeafOrder()
{
}

OptimalLeafOrder::~OptimalLeafOrder()
{
}

int 
OptimalLeafOrder::setHierarchicalClusterNode(HierarchicalClusterNode* aPtr)
{
	hcRoot=aPtr;
	return 0;
}

int
OptimalLeafOrder::reorder(vector<string>& ordering)
{
	root=new OptimalLeafOrder::Node;
	makeTree(hcRoot,root);
	if(root->left==NULL && root->right==NULL)
	{
		ordering.push_back(root->nodeName);
		return 0;
	}
	populateDist(hcRoot);
	reorder(root); 
	//Now get the highest scoring left and right child
	double min=10E7;
	string leftmost;
	string rightmost;
	for(map<string,map<string,double>*>::iterator cIter=root->childOrder.begin();cIter!=root->childOrder.end();cIter++)
	{
		map<string,double>* vals=cIter->second;
		for(map<string,double>::iterator vIter=vals->begin();vIter!=vals->end();vIter++)
		{
			//cout << cIter->first <<"\t" << vIter->first <<"\t" <<vIter->second << endl;
			if(vIter->second<min)
			{		
				min=vIter->second;
				leftmost.clear();
				rightmost.clear();
				leftmost.append(cIter->first);
				rightmost.append(vIter->first);
			}
		}
	}
	cout<< "Mindist " <<min << endl;
	string key;
	if(strcmp(leftmost.c_str(),rightmost.c_str())<=0)
	{
		key.append(leftmost.c_str());
		key.append("-");
		key.append(rightmost.c_str());	
	}
	else
	{
		key.append(rightmost.c_str());	
		key.append("-");
		key.append(leftmost.c_str());
	}
	//ordering.push_back(leftmost);
	OptimalLeafOrder::Pair* p=root->extremes[key];
	//getordering(root->left,leftmost,p->leftextreme,ordering);
	//getordering(root->right,p->rightextreme,rightmost,ordering);
	
	if(root->left->childOrder.find(leftmost)!=root->left->childOrder.end())	
	{
		getordering(root->left,leftmost,p->leftextreme,ordering);
		getordering(root->right,p->rightextreme,rightmost,ordering);
	}
	//switch case
	else if(root->right->childOrder.find(leftmost)!=root->right->childOrder.end())	
	{
		getordering(root->right,leftmost,p->rightextreme,ordering);
		getordering(root->left,p->leftextreme,rightmost,ordering);
	}
	double d=0;
	for(int i=0;i<ordering.size()-1;i++)
	{
		d=d+getSim(ordering[i],ordering[i+1]);
	}
	cout << "distance: " << d << endl;
	return 0;
}


int
OptimalLeafOrder::populateDist(HierarchicalClusterNode* node)
{
	if(node->left==NULL && node->right==NULL)
	{
		for(map<string,double>::iterator nIter=node->distToNeighbors.begin();nIter!=node->distToNeighbors.end();nIter++)
		{
			string key;
			if(strcmp(node->nodeName.c_str(),nIter->first.c_str())<=0)
			{
				key.append(node->nodeName.c_str());
				key.append("-");
				key.append(nIter->first.c_str());
			}
			else
			{
				key.append(nIter->first.c_str());
				key.append("-");
				key.append(node->nodeName.c_str());
			}
			if(nIter->second>10E7)
			{
				cout <<"Distance overflow error for "<<  key << " "<< nIter->second <<  endl;
				exit(0);
			}
			distance[key]=nIter->second;
		}
	}
	if(node->left!=NULL)
	{
		populateDist(node->left);
	}
	if(node->right!=NULL)
	{
		populateDist(node->right);
	}
	
	return 0;
}

int
OptimalLeafOrder::makeTree(HierarchicalClusterNode* node,OptimalLeafOrder::Node* ooNode)
{
	ooNode->left=NULL;
	ooNode->right=NULL;
	ooNode->parent=NULL;
	ooNode->nodeName.append(node->nodeName);
	if(node->left!=NULL)
	{
		OptimalLeafOrder::Node* leftchild=new OptimalLeafOrder::Node;
		ooNode->left=leftchild;
		makeTree(node->left,leftchild);
	}
	if(node->right!=NULL)
	{
		OptimalLeafOrder::Node* rightchild=new OptimalLeafOrder::Node;
		ooNode->right=rightchild;
		makeTree(node->right,rightchild);
	}
	return 0;
}

int
OptimalLeafOrder::reorder(OptimalLeafOrder::Node* n)
{
	//If this is left node then this is the only node possible here
	if(n->left==NULL && n->right==NULL)
	{
		map<string,double>* v=new map<string,double>;
		n->childOrder[n->nodeName]=v;
		(*v)[n->nodeName]=0;
		return 0;
	}
	reorder(n->left);
	reorder(n->right);
	map<string,map<string,double>*>& leftleaves=n->left->childOrder;
	map<string,map<string,double>*>& rightleaves=n->right->childOrder;
	for(map<string,map<string,double>*>::iterator vlIter=leftleaves.begin();vlIter!=leftleaves.end();vlIter++)
	{
		//Possible ms in vl
		map<string,double>* mset=vlIter->second; 
		for(map<string,map<string,double>*>::iterator ulIter=rightleaves.begin();ulIter!=rightleaves.end();ulIter++)
		{
			//Possible ks in vr
			map<string,double>* nset=ulIter->second;
			double mindist=10E7;
			string xleft;
			string xright;
			for(map<string,double>::iterator mIter=mset->begin();mIter!=mset->end();mIter++)
			{
				for(map<string,double>::iterator nIter=nset->begin();nIter!=nset->end();nIter++)
				{
					double score=getSim((string&)mIter->first,(string&)nIter->first)+mIter->second+nIter->second;
					if(score<mindist)
					{
						mindist=score;
						if(strcmp(mIter->first.c_str(),vlIter->first.c_str())!=0)
						{
							xleft.clear();
							xleft.append(mIter->first);
						}
						if(strcmp(nIter->first.c_str(),ulIter->first.c_str())!=0)
						{
							xright.clear();
							xright.append(nIter->first);
						}
					}
				}
			}
			if(mindist==1000)
			{
				cout <<"odd.. very odd. mindist is still small" <<endl;
			}
			map<string,double>* vals=NULL;
			if(n->childOrder.find(vlIter->first)==n->childOrder.end())
			{
				vals=new map<string,double>;
				n->childOrder[vlIter->first]=vals;
			}
			else
			{
				vals=n->childOrder[vlIter->first];
			}
			(*vals)[ulIter->first]=mindist;
			vals=NULL;
			if(n->childOrder.find(ulIter->first)==n->childOrder.end())
			{
				vals=new map<string,double>;
				n->childOrder[ulIter->first]=vals;
			}
			else
			{
				vals=n->childOrder[ulIter->first];
			}
			(*vals)[vlIter->first]=mindist;
			string key;
			if(strcmp(ulIter->first.c_str(),vlIter->first.c_str())<=0)
			{
				key.append(ulIter->first.c_str());
				key.append("-");
				key.append(vlIter->first.c_str());
			}
			else
			{
				key.append(vlIter->first.c_str());
				key.append("-");
				key.append(ulIter->first.c_str());
			}
			if(xleft.length()==0 && xright.length()==0)
			{
				continue;
			}
			OptimalLeafOrder::Pair* p=new OptimalLeafOrder::Pair;
			if(xleft.length()>0)
			{
				p->leftextreme.append(xleft.c_str());
			}
			if(xright.length()>0)
			{
				p->rightextreme.append(xright.c_str());
			}
			n->extremes[key]=p;
		}
	}
	return 0;
}

int
OptimalLeafOrder::getordering(OptimalLeafOrder::Node* n, string& leftextreme, string& rightextreme, vector<string>& leafordering)
{
	if(n->left==NULL && n->right==NULL)
	{
		leafordering.push_back(n->nodeName);
		return 0;
	}
	else if(n->extremes.size()==0)
	{
		leafordering.push_back(leftextreme);
		leafordering.push_back(rightextreme);
	}
	else
	{
		//leafordering.push_back(leftextreme);
		//leafordering.push_back(leftextreme);
		string key;
		if(strcmp(leftextreme.c_str(),rightextreme.c_str())<=0)
		{
			key.append(leftextreme.c_str());
			key.append("-");
			key.append(rightextreme.c_str());
		}
		else
		{
			key.append(rightextreme.c_str());
			key.append("-");
			key.append(leftextreme.c_str());
		}
		if(n->extremes.find(key)==n->extremes.end())
		{
			cout <<"No key " << key << " at " << n->nodeName <<  endl;
			exit(0);
		}

		Pair* ex=n->extremes[key];
		//if(ex->leftextreme.size()>0)
		//{
			//because the left and rights can switch, it is important to check that we are
			//setting the boundaries of the right trees
			if(n->left->childOrder.find(leftextreme)!=n->left->childOrder.end())	
			{
				getordering(n->left,leftextreme,ex->leftextreme,leafordering);
				getordering(n->right,ex->rightextreme,rightextreme,leafordering);
			}
			//switch case
			else if(n->right->childOrder.find(leftextreme)!=n->right->childOrder.end())	
			{
				getordering(n->right,leftextreme,ex->rightextreme,leafordering);
				getordering(n->left,ex->leftextreme,rightextreme,leafordering);
			}
		//}
		//if(ex->rightextreme.size()>0)
		//{
		//}
	}

	return 0;
}

double
OptimalLeafOrder::getSim(string& s1,string& s2)
{
	string key;
	if(strcmp(s1.c_str(),s2.c_str())<=0)
	{
		key.append(s1.c_str());
		key.append("-");
		key.append(s2.c_str());
	}
	else
	{
		key.append(s2.c_str());
		key.append("-");
		key.append(s1.c_str());
	}
	if(distance.find(key)==distance.end())
	{
		cout <<"No distance for " << key <<endl;
	}
	double d=distance[key];
	return d;
}
