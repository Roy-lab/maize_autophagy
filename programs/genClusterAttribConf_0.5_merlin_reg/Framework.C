#include <fstream>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include "Heap.H"
#include "HierarchicalClusterNode.H"
#include "HierarchicalCluster.H"
#include "OptimalLeafOrder.H"
#include "GeneNameMapper.H"
#include "GeneExpManager.H"
#include "Framework.H"
#include <unistd.h>

/*
 * DC updated this version to take an optional argument containing
 * confidence values for expression reg->gene edges.
 */
Framework::Framework()
{
	minT=0;
	maxT=1500;
}

Framework::~Framework()
{
}

int
Framework::init(int argc, char** argv)
{
	gnm.readGeneNames();
	int optret='-';
	opterr=1;
	int oldoptind=optind;
	int condCnt=1;
	char configFName[1024];	
	char chipEnrFName[1024];	
	char goEnrFName[1024];
	char expregEnrFName[1024];
	char dnaseEnrFName[1024];
	char confFName[1024];

	// initialize
	hasReadEdgeConfs=false;
	
	while(optret=getopt(argc,argv,"l:o:C:c:r:g:d:f:h:e:m:t:T:")!=-1)
	{
		if(optret=='?')
		{
			cout <<"Option error " << optopt << endl;
			return -1;
		}
		char c;
		char* my_optarg=NULL;
		c=*(argv[oldoptind]+1);
		if(optind-oldoptind ==2)
		{
			my_optarg=argv[oldoptind+1];	
		}
		else
		{
			my_optarg=argv[oldoptind]+2;
		}
		switch(c)
		{
			case 't':
			{
				minT = atoi(my_optarg);
				break;
			}
			case 'T':
			{
				maxT = atoi(my_optarg);
				break;
			}
			case 'l':
			{
				readSelectClusterEnrPairs(my_optarg);
				break;
			}
			case 'o':
			{
				strcpy(outDirName,my_optarg);
				break;
			}
			case 'C':
			{
				strcpy(configFName,my_optarg);
				break;
			}
			case 'c':
			{
				strcpy(chipEnrFName,my_optarg);
				break;
			}
			case 'r':
			{
				strcpy(expregEnrFName,my_optarg);
				break;
			}
			case 'g':
			{
				strcpy(goEnrFName,my_optarg);
				break;
			}
			case 'd':
			{
				strcpy(dnaseEnrFName,my_optarg);
				break;
			}
			case 'f':	 // DC ADDED
			{
				strcpy(confFName,my_optarg);
				break;
			}
			case 'h':
			{
				readHeaderNames(my_optarg);
				break;
			}
			case 'e':
			{
				expMgr.readExpression(my_optarg);
				break;
			}
			case 'm':
			{
				readClusterMembership(my_optarg);
				break;
			}
			default:
			{
				cout <<"Unhandled option " << c  << endl;
				return -1;
			}
		}
		oldoptind=optind;
	}

	int success=0;
	if(strlen(configFName)>0)
	{
		readConfig(configFName);
	}
	else
	{
		if(strlen(chipEnrFName)>0)
		{
			string enrType("3");
			readEnr(chipEnrFName,enrType);
		}
		if(strlen(expregEnrFName)>0)
		{
			string enrType("2");
			readEnr(expregEnrFName,enrType);
		}
		if(strlen(dnaseEnrFName)>0)
		{
			string enrType("4");
			readEnr(dnaseEnrFName,enrType);
		}
		if(strlen(goEnrFName)>0)
		{
			string enrType("1");
			readEnr(goEnrFName,enrType);
		}
	}
	// DC ADDED -- read edge confs
	if (strlen(confFName)>0)
	{
		success=readEdgeConfs(confFName);
		hasReadEdgeConfs=true;
	}
	//filterClusters();

	return success;
}

/*
 * Generate the _attrib.txt file
 * Updated by DC to optionally use edge confs instead of 0/1 for enriched regs.
 */
int
Framework::genClusterAttribs()
{
	map<string, vector<double>* > modAvgMap;
	for(map<int,map<string,int>*>::iterator cIter=clusterEnrPair.begin();cIter!=clusterEnrPair.end();cIter++)
	{
		map<string,int>* cmembers=exprClust[cIter->first];

		if (cmembers->size() < minT || cmembers->size()>maxT)
		{
			cout << "Skipping cluster " << cIter->first << " of size " << cmembers->size() << endl;
			continue;
		}

		cout <<"Analyzing cluster "<< cIter->first << " with " << cmembers->size() << " members"<< endl;
		if(cIter->first==0)
		{
			cout <<"Stop here" << endl;
		}
		char outFName[1024];
		sprintf(outFName,"%s/Cluster%d_attrib.txt",outDirName,cIter->first);
		map<string,int>* attrib=cIter->second;
		map<int,map<string,int>*> grouping;
		for(map<string,int>::iterator aIter=attrib->begin();aIter!=attrib->end();aIter++)
		{
			
			char termName[1024];
			strcpy(termName,aIter->first.c_str());
			char* pos=strrchr(termName,':');
			int type;
			if(pos!=NULL)
			{
				*pos='\0';
				type=atoi(pos+1);
			}
			else
			{	
				cout <<"Bad format no colon!" << endl;
				exit(0);
			}
			map<string,int>* grp=NULL;
			if(grouping.find(type)==grouping.end())
			{
				grp=new map<string,int>;
				grouping[type]=grp;
			}
			else
			{
				grp=grouping[type];
			}
			(*grp)[termName]=type;
		}
		vector<string> colAttribs;
		vector<string> rowAttribs;
		//For now I am going to ignore a module that does not have a expression regulators
		//if(grouping.find(2)==grouping.end())
		//{
		//	cout <<"No expression regulators for cluster " << cIter->first << endl;
		//	continue;
		//}
		//Order by expression regulators -- grouping value 2
		if(grouping.find(2)!=grouping.end())
		{
			if(grouping[2]->size()>2)
			{
				orderColumns_Optimal(cIter->first,grouping[2],2,colAttribs);
			}
			else
			{
				orderColumns(cIter->first,grouping[2],2,colAttribs);
			}
			if(cmembers->size()>200)
			{
				orderGenes(cIter->first,grouping[2],2,rowAttribs);
				//orderGenes_Optimal(cIter->first,grouping[2],2,rowAttribs);
			}
			else
			{
				orderGenes_Optimal(cIter->first,grouping[2],2,rowAttribs);
			}
		}	
		else if(grouping.find(3)!=grouping.end()) // ChIP
		{
			if(grouping[3]->size()>2)
			{
				orderColumns_Optimal(cIter->first,grouping[3],3,colAttribs);
			}
			else
			{
				orderColumns(cIter->first,grouping[3],3,colAttribs);
			}
			orderGenes_Optimal(cIter->first,grouping[3],3,rowAttribs);
		}	
		else
		{	
			orderGenes_Optimal(cIter->first, grouping[6], 6, rowAttribs);
			cout <<"No ChIP or Expression regulators in cluster " << cIter->first << endl;
			
			//continue;
		}
		vector<string> gocolAttribs;
		if(grouping.find(1)!=grouping.end())
		{
			if(grouping[1]->size()>2)
			{
				orderColumns_Optimal(cIter->first,grouping[1],1,gocolAttribs);
			}
			else
			{
				orderColumns(cIter->first,grouping[1],1,gocolAttribs);
			}
		}
		/*else
		{
			int gid=grouping.begin()->first;
			cout <<"No Go process ordering with "<< gid << endl;
			orderColumns(cIter->first,grouping[gid],gid,gocolAttribs);
			//continue;
		}*/
		ofstream oFile(outFName);
	
		//orderGenes(cIter->first,grouping[2],2,rowAttribs);
		vector< vector<double>* > moduleExp;
		for(int g=0;g<rowAttribs.size();g++)
		{
			cout << g << "\t" << rowAttribs[g] << endl;
			//first show expression
			//Now the other attributes
			map<string,string>* geneAttrib=NULL;
			if(geneEnrMap.find(rowAttribs[g])!=geneEnrMap.end())	
			{
				geneAttrib=geneEnrMap[rowAttribs[g]];
			}
			for(map<int,map<string,int>*>::iterator aIter=grouping.begin();aIter!=grouping.end();aIter++)
			{
				if(aIter->first==1)
				{
					for(int j=0;j<gocolAttribs.size();j++)
					{
						char tempName[1024];
						sprintf(tempName,"%s:%d",gocolAttribs[j].c_str(),aIter->first);

						string key;
						key.append(tempName);
						string value("0");
						if(geneAttrib!=NULL)
						{
							if(geneAttrib->find(key)!=geneAttrib->end())
							{
								value.clear();
								value.append((*geneAttrib)[key]);
							}
						}
						oFile << gnm.getCommonName(rowAttribs[g].c_str()) << "||6\t" << key << "||6\t" << value<<"|1" << endl;
					}
				}
				else if(aIter->first==3 || aIter->first==4 || aIter->first==5)
				{
					map<string,int>* attribname=aIter->second;
					for(map<string,int>::iterator tIter=attribname->begin();tIter!=attribname->end();tIter++)
					{
						char tempName[1024];
						sprintf(tempName,"%s:%d",tIter->first.c_str(),aIter->first);
						string key;
						key.append(tempName);
						string value("0");
						if(geneAttrib!=NULL)
						{
							if(geneAttrib->find(key)!=geneAttrib->end())
							{
								value.clear();
								value.append((*geneAttrib)[key]);
							}
						}
						if((aIter->first==2) || (aIter->first==3) || (aIter->first==4)|| (aIter->first==5))
						{
							oFile << gnm.getCommonName(rowAttribs[g].c_str()) << "||6\t" << gnm.getCommonName(tIter->first.c_str()) << ":"<< aIter->first<< "||6\t" << value<<"|1" << endl;
						}
						else
						{
							oFile << gnm.getCommonName(rowAttribs[g].c_str()) << "||6\t" << key << "||6\t" << value<<"|1" << endl;
						}
					}
				}
				else if(aIter->first==2) // expression regulators
				{
					for(int j=0;j<colAttribs.size();j++)
					{
						char tempName[1024];
						sprintf(tempName,"%s:%d",colAttribs[j].c_str(),aIter->first);

						string key;
						key.append(tempName);
						string value("0");
						if(geneAttrib!=NULL)
						{
							if(geneAttrib->find(key)!=geneAttrib->end())
							{
								value.clear();
								value.append((*geneAttrib)[key]);
							}

							if(hasReadEdgeConfs) // Get conf val for this edge
							{
								value.clear();
								string myKey; // reg-target
								myKey.append(colAttribs[j].c_str());
								myKey.append("-");
								myKey.append(rowAttribs[g].c_str());
								char tempVal[1024];
								sprintf(tempVal, "%f", edgeConfs[myKey]);
								value.append(tempVal);
							}

						}
						oFile << gnm.getCommonName(rowAttribs[g].c_str()) << "||6\t" << gnm.getCommonName(colAttribs[j].c_str()) << ":"<< aIter->first<< "||6\t" << value<<"|1" << endl;
					}
				}
				if(aIter->first != 6)
				{
					if(g==0)
					{
						oFile <<"|- "<<"VertSpacer"<< aIter->first <<"|Spacer|3" << endl;
					}
				}
			}
			// now show expression
			//cout << "showing " << rowAttribs[g] << endl;
			vector<double>* expr=expMgr.getExp(rowAttribs[g]);
			for(int i=0;i<expr->size();i++)
			{
				oFile << gnm.getCommonName(rowAttribs[g].c_str())  <<"||6\t" << headerNames[i]<<"||6\t" << (*expr)[i] <<"|2" << endl;
			}
			moduleExp.push_back(expr);
		}
		if (moduleExp.size()>1)
		{
			//oFile <<"|Spacer|3|HorSpaceModuleExp|-"<<  endl;	
			vector<double>* modAvg = getAverage(moduleExp);
			//for(int i=0;i<modAvg->size();i++)
			//{
			//	oFile << "moduleAvg"  <<"||6\t" << headerNames[i]<<"||6\t" << (*modAvg)[i] <<"|2" << endl;
			//}
			char cName[1024];
			sprintf(cName,"Cluster%d",cIter->first);
			string cNameStr=string(cName);
			modAvgMap[cNameStr] = modAvg;
			//delete modAvg;
		}
		moduleExp.clear();

		oFile.close();
	}
	char modAvgFName[1024];
	vector<string> modOrder;
	sprintf(modAvgFName,"%s/ModuleAvg.txt",outDirName);
	ofstream oModFile(modAvgFName);
	orderModAvg_Optimal(modAvgMap, modOrder);
	for (int j=0;j<modOrder.size();j++)
	{
		string cName=modOrder[j];
		vector<double>* expr=modAvgMap[cName];
		for(int i=0;i<expr->size();i++)
		{
			//oModFile << cName <<"||6\t" << headerNames[i]<<"||6\t" << (*expr)[i] <<"|2" << endl;
			oModFile << cName <<"||6\t" << headerNames[i]<<"||6\t" << (*expr)[i] <<"|" << endl;
		}
	}
	oModFile.close();
	return 0;
}

vector<double>*
Framework::getAverage(vector< vector<double>* >& moduleExp)
{
	vector<double>* res = new vector<double>;
	for (int j=0;j<moduleExp[0]->size();j++)
	{
		res->push_back(moduleExp[0]->at(j));
	}
	for (int i=1;i<moduleExp.size();i++)
	{
		for (int j=0;j<res->size();j++)
		{
			(*res)[j] = res->at(j) + moduleExp[i]->at(j);
		}
	}
	int modSize = double(moduleExp.size());
	for (int j=0;j<res->size();j++)
	{
		(*res)[j] = res->at(j)/modSize;
	}
	return res;
}

/*int
Framework::genClusterAttribs()
{
	for(map<int,map<string,int>*>::iterator cIter=clusterEnrPair.begin();cIter!=clusterEnrPair.end();cIter++)
	{
		cout <<"Analyzing cluster "<< cIter->first << endl;
		if(cIter->first==107)
		{
			cout <<"Stop here" << endl;
		}
		char outFName[1024];
		sprintf(outFName,"%s/Cluster%d_attrib.txt",outDirName,cIter->first);
		ofstream oFile(outFName);
		map<string,int>* attrib=cIter->second;
		map<string,int>* cmembers=exprClust[cIter->first];
		map<int,map<string,int>*> grouping;
		for(map<string,int>::iterator aIter=attrib->begin();aIter!=attrib->end();aIter++)
		{
			
			char termName[1024];
			strcpy(termName,aIter->first.c_str());
			char* pos=strrchr(termName,':');
			int type;
			if(pos!=NULL)
			{
				*pos='\0';
				type=atoi(pos+1);
			}
			else
			{	
				cout <<"Bad format no colon!" << endl;
				exit(0);
			}
			map<string,int>* grp=NULL;
			if(grouping.find(type)==grouping.end())
			{
				grp=new map<string,int>;
				grouping[type]=grp;
			}
			else
			{
				grp=grouping[type];
			}
			(*grp)[termName]=type;
		}
		vector<string> colAttribs;
		vector<string> rowAttribs;
		//Order by expression regulators 
		orderColumns(cIter->first,grouping[2],2,colAttribs);
		orderGenes_Cnt(cIter->first,grouping[2],2,rowAttribs);
		for(map<int,map<string,int>*>::iterator aIter=grouping.begin();aIter!=grouping.end();aIter++)
		{
			for(int g=0;g<rowAttribs.size();g++)
			{
				//first show expression
				//Now the other attributes
				map<string,string>* geneAttrib=NULL;
				if(geneEnrMap.find(rowAttribs[g])!=geneEnrMap.end())	
				{
					geneAttrib=geneEnrMap[rowAttribs[g]];
				}
				if(aIter->first!=2)
				{
					map<string,int>* attribname=aIter->second;
					for(map<string,int>::iterator tIter=attribname->begin();tIter!=attribname->end();tIter++)
					{
						char tempName[1024];
						sprintf(tempName,"%s:%d",tIter->first.c_str(),aIter->first);
						string key;
						key.append(tempName);
						string value("0");
						if(geneAttrib!=NULL)
						{
							if(geneAttrib->find(key)!=geneAttrib->end())
							{
								value.clear();
								value.append((*geneAttrib)[key]);
							}
						}
						if((aIter->first==2) || (aIter->first==3))
						{
							oFile << gnm.getCommonName(rowAttribs[g].c_str()) << "||6\t" << gnm.getCommonName(tIter->first.c_str()) << ":"<< aIter->first<< "||6\t" << value<<"|1" << endl;
						}
						else
						{
							oFile << gnm.getCommonName(rowAttribs[g].c_str()) << "||6\t" << key << "||6\t" << value<<"|1" << endl;
						}
					}
				}
				else
				{
					for(int j=0;j<colAttribs.size();j++)
					{
						char tempName[1024];
						sprintf(tempName,"%s:%d",colAttribs[j].c_str(),aIter->first);

						string key;
						key.append(tempName);
						string value("0");
						if(geneAttrib!=NULL)
						{
							if(geneAttrib->find(key)!=geneAttrib->end())
							{
								value.clear();
								value.append((*geneAttrib)[key]);
							}
						}
						oFile << gnm.getCommonName(rowAttribs[g].c_str()) << "||6\t" << gnm.getCommonName(colAttribs[j].c_str()) << ":"<< aIter->first<< "||6\t" << value<<"|1" << endl;
					}
				}
			}
		}
		oFile <<"|- VertSpacer|Spacer" << endl;
		for(int g=0;g<rowAttribs.size();g++)
		{
			vector<double>* expr=expMgr.getExp(rowAttribs[g]);
			for(int i=0;i<expr->size();i++)
			{
				oFile << gnm.getCommonName(rowAttribs[g].c_str())  <<"||6\t" << headerNames[i]<<"||6\t" << (*expr)[i] <<"|2" << endl;
			}
		}

		oFile.close();
	}
	return 0;
}*/



int
Framework::genClusterRegulators()
{
	for(map<int,map<string,int>*>::iterator cIter=clusterEnrPair.begin();cIter!=clusterEnrPair.end();cIter++)
	{
		map<string,int>* cmembers = exprClust[cIter->first];
		if (cmembers->size() < minT || cmembers->size()>maxT)
		{
			cout << "Skipping cluster " << cIter->first << " of size " << cmembers->size() << endl;
			continue;
		}

		char outFName[1024];
		sprintf(outFName,"%s/Cluster%d_regulators.txt",outDirName,cIter->first);
		ofstream oFile(outFName);
		map<string,int>* attrib=cIter->second;
		map<int,map<string,int>*> grouping;
		for(map<string,int>::iterator aIter=attrib->begin();aIter!=attrib->end();aIter++)
		{
			
			char termName[1024];
			strcpy(termName,aIter->first.c_str());
			char* pos=strrchr(termName,':');
			int type;
			if(pos!=NULL)
			{
				*pos='\0';
				type=atoi(pos+1);
			}
			else
			{	
				cout <<"Bad format no colon!" << endl;
				exit(0);
			}
			map<string,int>* grp=NULL;
			if(grouping.find(type)==grouping.end())
			{
				grp=new map<string,int>;
				grouping[type]=grp;
			}
			else
			{
				grp=grouping[type];
			}
			(*grp)[termName]=type;
		}
		//Show Expression regulators
		if(grouping.find(2)==grouping.end())
		{
			cout <<"No expression based regulators! " << endl;
			continue;
		}
		map<string,int>* exprReg=grouping[2];
		oFile <<"|Spacer|3|HorSpaceE"<<cIter->first << "|-"<<  endl;	
		for(map<string,int>::iterator rIter=exprReg->begin();rIter!=exprReg->end();rIter++)
		{
			vector<double>* expr=expMgr.getExp(rIter->first);
			if(expr==NULL)
			{
				continue;
			}
			for(int i=0;i<expr->size();i++)
			{
				oFile << gnm.getCommonName(rIter->first.c_str()) <<":2" <<"||6\t" << headerNames[i]<<"||6\t" << (*expr)[i] <<"|2" << endl;
			}
		}
		if(grouping.find(3)==grouping.end())
		{
			cout <<"No chip based regulators! " << endl;
		}
		else
		{
			map<string,int>* chipReg=grouping[3];
			oFile <<"|Spacer|3|HorSpaceC"<<cIter->first << "|-"<<  endl;	
			for(map<string,int>::iterator rIter=chipReg->begin();rIter!=chipReg->end();rIter++)
			{
				vector<double>* expr=expMgr.getExp(rIter->first);
				if(expr==NULL)
				{
					continue;
				}
				for(int i=0;i<expr->size();i++)
				{
					oFile << gnm.getCommonName(rIter->first.c_str()) <<":3" <<"||6\t" << headerNames[i]<<"||6\t" << (*expr)[i] <<"|2" << endl;
				}
			}
		}
		if(grouping.find(4)==grouping.end())
		{
			cout <<"No motif based regulators! " << endl;
			continue;
		}
		map<string,int>* motifReg=grouping[4];
		oFile <<"|Spacer|3|HorSpaceM"<<cIter->first << "|-"<<  endl;	
		for(map<string,int>::iterator rIter=motifReg->begin();rIter!=motifReg->end();rIter++)
		{
			vector<double>* expr=expMgr.getExp(rIter->first);
			if(expr==NULL)
			{
				continue;
			}
			for(int i=0;i<expr->size();i++)
			{
				oFile << gnm.getCommonName(rIter->first.c_str()) <<":4" <<"||6\t" << headerNames[i]<<"||6\t" << (*expr)[i] <<"|2" << endl;
			}
		}
			
		for(map<int,map<string,int>*>::iterator aIter=grouping.begin();aIter!=grouping.end();aIter++)
		{
			map<string,int>* attribname=aIter->second;
			attribname->clear();
			delete attribname;
		}
		grouping.clear();
		oFile.close();
	}
	return 0;
}

int 
Framework::readSelectClusterEnrPairs(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1024);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		int tokCnt=0;
		int cid=0;
		string enrType;
		string enrName;
		char* tok=strtok(buffer,"\t");
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				enrName.append(tok);
			}
			else if(tokCnt==1)
			{
				cid=atoi(tok+strlen("Cluster"));
			}
			else if(tokCnt==2)
			{
				enrType.append(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		map<string,int>* enrSet=NULL;
		if(clusterEnrPair.find(cid)==clusterEnrPair.end())
		{
			enrSet=new map<string,int>;
			clusterEnrPair[cid]=enrSet;
		}
		else 
		{
			enrSet=clusterEnrPair[cid];
		}
		string key(enrName);
		key.append(":");
		key.append(enrType);
		(*enrSet)[key]=0;
	}
	inFile.close();
	return 0;
}

int 
Framework::readEnr(const char* aFName, string& enrtype)
{
	ifstream inFile(aFName);
	char* buffer=NULL;
	char* memberSet=NULL;
	string buffstr;
	int bufflen=0;
	while(inFile.good())
	{
		getline(inFile,buffstr);	
		if(buffstr.length()<=0)
		{
			continue;
		}
		if(bufflen<=buffstr.length())
		{
			bufflen=buffstr.length()+1;
			if(buffer!=NULL)
			{
				delete[] buffer;
				delete[] memberSet;
			}
			buffer=new char [bufflen];
			memberSet=new char [bufflen];
		}
		strcpy(buffer,buffstr.c_str());
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		int cid=0;
		string term;	
		while(tok!=NULL)
		{
			switch(tokCnt)
			{
				case 0:
				{
					cid=atoi(tok+strlen("Cluster"));
					break;
				}
				case 1:
				{
					int i=0;
					while(tok[i]!='\0')
					{	
						if(isspace(tok[i]))
						{
							tok[i]='_';
						}
						i++;
					}
					term.append(tok);
					break;
				}
				case 9:
				{
					strcpy(memberSet,tok);
					break;
				}
				default:
				{
					break;
				}
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		if(clusterEnrPair.find(cid)==clusterEnrPair.end())
		{
			continue;
		}
		string key(term);
		key.append(":");
		key.append(enrtype);
		if(strcmp(term.c_str(),"YOR344C")==0 && cid==107) 
		{
			cout << "Stop here " << endl;
		}
		map<string,int>* enrSet=clusterEnrPair[cid];
		if(enrSet->find(key)==enrSet->end())
		{
			continue;
		}
		char* gene=strtok(memberSet,";");
		string geneName;
		while(gene!=NULL)
		{
			string orfName(gnm.getORFName(gene));
			map<string,string>* enrToShowForGene=NULL;
			if(geneEnrMap.find(orfName)==geneEnrMap.end())
			{
				enrToShowForGene=new map<string,string>;
				geneEnrMap[orfName]=enrToShowForGene;
			}
			else
			{
				enrToShowForGene=geneEnrMap[orfName];
			}
			(*enrToShowForGene)[key]=enrtype;
			gene=strtok(NULL,";");
		}
	}
	inFile.close();
	return 0;
}

int 
Framework::readClusterMembership(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good())	
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string geneName;
		int cassign;
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				geneName.append(tok);
			}
			else
			{
				cassign=atoi(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		map<string,int>* clustMembers=NULL;
		if(exprClust.find(cassign)==exprClust.end())
		{
			clustMembers=new map<string,int>;
			exprClust[cassign]=clustMembers;
		}
		else
		{
			clustMembers=exprClust[cassign];
		}
		(*clustMembers)[geneName]=0;
	}
	inFile.close();
	return 0;
}

int 
Framework::readHeaderNames(const char* headerFName)
{
	ifstream inFile(headerFName);
	int id=0;
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string headerName;
		headerName.append(buffer);
		headerNames[id]=headerName;
		id++;
	}
	inFile.close();
	return 0;
}

/*
 *DC adding -- read in edge confs to visualize instead of 0/1 for inferred regs.
 * We took enrichment from a low conf value, but will visualize the conf value for the edge here.
 * Input format in three cols: reg, target, conf value.
 */
int 
Framework::readEdgeConfs(const char* confFName)
{
	ifstream infile(confFName);
	char buffer[1024];
	while(infile.good())
	{
		infile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string key;
		double conf; 

		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				key.append(tok);
			}
			else if (tokCnt==1)
			{
				key.append("-");
				key.append(tok);
			} 
			else if (tokCnt==2)
			{
				conf=atof(tok);
			}	
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		//cout << key << ": " << conf << endl;
		if (edgeConfs.find(key) != edgeConfs.end() )
		{
			cerr << "Duplicate edge conf found: " << key << " " << conf << " " << edgeConfs[key] << endl;
			return 1;
		}
		edgeConfs[key]=conf;
	}
	infile.close();
	cout << "Read " << edgeConfs.size() << " edge conf values." << endl;
	return 0;	
}

int 
Framework::orderColumns(int cID,map<string,int>* attribNames, int attribKey, vector<string>& attribOrder)
{	
	map<string,int>* members=exprClust[cID];
	map<string,int> tfHits;
	for(map<string,int>::iterator aIter=attribNames->begin();aIter!=attribNames->end();aIter++)
	{
		char buffer[1024];
		sprintf(buffer,"%s:%d",aIter->first.c_str(), attribKey);	
		string key(buffer);
		for(map<string,int>::iterator mIter=members->begin();mIter!=members->end();mIter++)
		{
			if(geneEnrMap.find(mIter->first)==geneEnrMap.end())
			{
				continue;
			}
			map<string,string>* geneattrib=geneEnrMap[mIter->first];
			if(geneattrib->find(key)==geneattrib->end())
			{
				continue;
			}
			if(tfHits.find(aIter->first)==tfHits.end())
			{
				tfHits[aIter->first]=1;
			}
			else
			{
				tfHits[aIter->first]=tfHits[aIter->first]+1;
			}
		}
	}
	for(map<string,int>::iterator tIter=tfHits.begin();tIter!=tfHits.end();tIter++)
	{
		attribOrder.push_back(tIter->first);		
	}
	for(int i=0;i<attribOrder.size();i++)
	{
		for(int j=i+1;j<attribOrder.size();j++)
		{
			int u=tfHits[attribOrder[i]];
			int v=tfHits[attribOrder[j]];
			if(u<v)
			{
				string temp(attribOrder[i]);
				attribOrder[i].clear();
				attribOrder[i].append(attribOrder[j].c_str());
				attribOrder[j].clear();
				attribOrder[j].append(temp.c_str());	
			}	
		}
	}
	return 0;
}

int
Framework::orderColumns_Optimal(int cID,map<string,int>* attribNames, int attribKey, vector<string>& attribOrder)
{
	map<string,int>* members=exprClust[cID];
	HierarchicalCluster c;
	map<string,HierarchicalClusterNode*> nodeSet;
	map<string,HierarchicalClusterNode*> backup;
	OptimalLeafOrder olo;
	map<string,int> attribNameIDMap;
	for(map<string,int>::iterator aIter=attribNames->begin();aIter!=attribNames->end();aIter++)
	{
		int attrcnt=0;
		HierarchicalClusterNode* hcNode=new HierarchicalClusterNode;
		nodeSet[aIter->first]=hcNode;
		backup[aIter->first]=hcNode;
		hcNode->nodeName.append(aIter->first.c_str());
		char buffer[1024];
		sprintf(buffer,"%s:%d",aIter->first.c_str(), attribKey);	
		string key(buffer);
		for(map<string,int>::iterator mIter=members->begin();mIter!=members->end();mIter++)
		{
			if(geneEnrMap.find(mIter->first)==geneEnrMap.end())
			{
				continue;
			}
			map<string,string>* geneattrib=geneEnrMap[mIter->first];
			if(geneattrib->find(key)==geneattrib->end())
			{
				continue;
			}
			int attribID=-1;
			if(attribNameIDMap.find(mIter->first)==attribNameIDMap.end())
			{
				attribID=attribNameIDMap.size();
				attribNameIDMap[mIter->first]=attribID;
			}
			else
			{
				attribID=attribNameIDMap[mIter->first];
			}
			double val=atof((*geneattrib)[key].c_str());
			if(val>0)
			{
				val=1;
			}
			hcNode->attrib[attribID]=val;
		}
	}
	map<int,map<string,int>*> modules;
	c.cluster(modules,nodeSet,0);
	olo.setHierarchicalClusterNode(c.getRoot());
	olo.reorder(attribOrder);
	return 0;
}


int 
Framework::orderGenes(int cID,map<string,int>* attribNames,int attribKey, vector<string>& geneOrder)
{
	map<string,int>* members=exprClust[cID];
	map<string,int> geneAttribs;
	string geneWithMostAttribs;
	int maxattrcnt=0;
	for(map<string,int>::iterator mIter=members->begin();mIter!=members->end();mIter++)
	{
		int attrcnt=0;
		for(map<string,int>::iterator aIter=attribNames->begin();aIter!=attribNames->end();aIter++)
		{
			char buffer[1024];
			sprintf(buffer,"%s:%d",aIter->first.c_str(), attribKey);	
			string key(buffer);
			if(geneEnrMap.find(mIter->first)==geneEnrMap.end())
			{
				continue;
			}
			map<string,string>* geneattrib=geneEnrMap[mIter->first];
			if(geneattrib->find(key)==geneattrib->end())
			{
				continue;
			}
			attrcnt++;
		}
		geneAttribs[mIter->first]=attrcnt;
		if(attrcnt>maxattrcnt)
		{
			geneWithMostAttribs.clear();
			geneWithMostAttribs.append(mIter->first.c_str());
			maxattrcnt=attrcnt;
		}
	}
	map<string,double> geneDist;
	geneOrder.push_back(geneWithMostAttribs.c_str());
	geneDist[geneWithMostAttribs]=1;
	map<string,string>* attrib1=geneEnrMap[geneWithMostAttribs];
	//Now compute similarity wrt to the genewithMaxAttr and order everything wrt it
	for(map<string,int>::iterator mIter=members->begin();mIter!=members->end();mIter++)
	{
		if(strcmp(mIter->first.c_str(),geneWithMostAttribs.c_str())==0)
		{
			continue;
		}
		int t1=0;
		int t2=0;
		int hit=0;
		if(geneEnrMap.find(mIter->first)==geneEnrMap.end())
		{
			geneDist[mIter->first]=0;
			geneOrder.push_back(mIter->first);
			continue;
		}
		map<string,string>* attrib=geneEnrMap[mIter->first];
		for(map<string,int>::iterator aIter=attribNames->begin();aIter!=attribNames->end();aIter++)
		{
			char keybuffer[1024];
			sprintf(keybuffer,"%s:%d",aIter->first.c_str(),attribKey);
			string key(keybuffer);
			if(attrib1->find(key)!=attrib1->end())
			{
				t1++;
			}
			if(attrib->find(key)!=attrib->end())
			{
				t2++;
			}
			if(attrib->find(key)!=attrib->end() && attrib1->find(key)!=attrib1->end())
			{
				hit++;
			}
		}
		double score=((double)hit)/((double)(t1+t2-hit));
		geneDist[mIter->first]=score;
		geneOrder.push_back(mIter->first);
	}
	for(int i=0;i<geneOrder.size();i++)
	{
		for(int j=i+1;j<geneOrder.size();j++)
		{
			double u=geneDist[geneOrder[i]];
			double v=geneDist[geneOrder[j]];
			if(u<v)
			{
				string temp(geneOrder[i]);
				geneOrder[i].clear();
				geneOrder[i].append(geneOrder[j].c_str());
				geneOrder[j].clear();
				geneOrder[j].append(temp.c_str());	
			}	
		}
	}
	return 0;
}


int 
Framework::orderGenes_Cnt(int cID,map<string,int>* attribNames,int attribKey, vector<string>& geneOrder)
{
	map<string,int>* members=exprClust[cID];
	map<string,int> geneAttribs;
	string geneWithMostAttribs;
	int maxattrcnt=0;
	for(map<string,int>::iterator mIter=members->begin();mIter!=members->end();mIter++)
	{
		int attrcnt=0;
		for(map<string,int>::iterator aIter=attribNames->begin();aIter!=attribNames->end();aIter++)
		{
			char buffer[1024];
			sprintf(buffer,"%s:%d",aIter->first.c_str(), attribKey);	
			string key(buffer);
			if(geneEnrMap.find(mIter->first)==geneEnrMap.end())
			{
				continue;
			}
			map<string,string>* geneattrib=geneEnrMap[mIter->first];
			if(geneattrib->find(key)==geneattrib->end())
			{
				continue;
			}
			attrcnt++;
		}
		geneAttribs[mIter->first]=attrcnt;
		if(attrcnt>maxattrcnt)
		{
			geneWithMostAttribs.clear();
			geneWithMostAttribs.append(mIter->first.c_str());
			maxattrcnt=attrcnt;
		}
	}
	//Now compute similarity wrt to the genewithMaxAttr and order everything wrt it
	for(map<string,int>::iterator mIter=members->begin();mIter!=members->end();mIter++)
	{
		geneOrder.push_back(mIter->first);
	}
	for(int i=0;i<geneOrder.size();i++)
	{
		for(int j=i+1;j<geneOrder.size();j++)
		{
			int u=geneAttribs[geneOrder[i]];
			int v=geneAttribs[geneOrder[j]];
			if(u<v)
			{
				string temp(geneOrder[i]);
				geneOrder[i].clear();
				geneOrder[i].append(geneOrder[j].c_str());
				geneOrder[j].clear();
				geneOrder[j].append(temp.c_str());	
			}	
		}
	}
	return 0;
}


int 
Framework::orderGenes_Optimal(int cID, map<string,int>* attribNames,int attribKey, vector<string>& geneOrder)
{
	map<string,int>* members=exprClust[cID];
	HierarchicalCluster c;
	map<string,HierarchicalClusterNode*> nodeSet;
	map<string,HierarchicalClusterNode*> backup;
	OptimalLeafOrder olo;
	map<string,int> attribNameIDMap;
	for(map<string,int>::iterator mIter=members->begin();mIter!=members->end();mIter++)
	{
		int attrcnt=0;
		HierarchicalClusterNode* hcNode=new HierarchicalClusterNode;
		nodeSet[mIter->first]=hcNode;
		backup[mIter->first]=hcNode;
		hcNode->nodeName.append(mIter->first.c_str());
		for(map<string,int>::iterator aIter=attribNames->begin();aIter!=attribNames->end();aIter++)
		{
			char buffer[1024];
			sprintf(buffer,"%s:%d",aIter->first.c_str(), attribKey);	
			string key(buffer);
			if(geneEnrMap.find(mIter->first)==geneEnrMap.end())
			{
				continue;
			}
			map<string,string>* geneattrib=geneEnrMap[mIter->first];
			if(geneattrib->find(key)==geneattrib->end())
			{
				continue;
			}
			int attribID=-1;
			if(attribNameIDMap.find(key)==attribNameIDMap.end())
			{
				attribID=attribNameIDMap.size();
				attribNameIDMap[key]=attribID;
			}
			else
			{
				attribID=attribNameIDMap[key];
			}
			double val=atof((*geneattrib)[key].c_str());
			if(val>0)
			{
				val=1;
			}
			hcNode->attrib[attribID]=val;
		}
	}
	map<int,map<string,int>*> modules;
	c.cluster(modules,nodeSet,0);
	olo.setHierarchicalClusterNode(c.getRoot());
	olo.reorder(geneOrder);
	return 0;
}


int 
Framework::orderModAvg_Optimal(map<string, vector<double>* >& modAvgMap, vector<string>& modOrder)
{
	HierarchicalCluster c;
	map<string,HierarchicalClusterNode*> nodeSet;
	map<string,HierarchicalClusterNode*> backup;
	OptimalLeafOrder olo;
	for(map<string, vector<double>* >::iterator mIter=modAvgMap.begin();mIter!=modAvgMap.end();mIter++)
	{
		vector<double>* expr = mIter->second;

		HierarchicalClusterNode* hcNode=new HierarchicalClusterNode;
		nodeSet[mIter->first]=hcNode;
		hcNode->nodeName = mIter->first;
		for (int i=0;i<expr->size();i++)
		{
			hcNode->expr.push_back(expr->at(i));
		}
	}
	map<int,map<string,int>*> modules;
	c.clusterExp(modules,nodeSet,0);
	olo.setHierarchicalClusterNode(c.getRoot());
	olo.reorder(modOrder);
	return 0;
}

/* Added by Ali
 * Instead of reading enrichment files one by one, we use 1 config file
 * It is: ID[\t]Name
 * Like:
 * 1	example/goproc_details.txt
*/
int
Framework::readConfig(const char* confFName)
{
	ifstream infile(confFName);
	char inname[1024];
	char buffer[1024];
	while(infile.good())
	{
		infile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string ID;

		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				ID=string(tok);
			}
			else if (tokCnt==1)
			{
				strcpy(inname,tok);
			} 
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		readEnr(inname,ID);
	}
	infile.close();
	return 0;
}

int
main(int argc, char* argv[])
{
	if(argc<2)
	{
		cout <<"genClusterAttrib" <<  endl
			<< "-l genelist"<< endl
			<< "-o outputdir" << endl
			<< "-t min cluster size (default 0)" << endl
			<< "-T max cluster size (default 1500)" << endl
			<< "-f conf values for exprreg edges" << endl
			<< "-h exprheader" << endl
			<< "-e exprmat" << endl
			<< "-m modulemembership" << endl
			<< "-C configFile" << endl
			<< "or alternatively:" << endl
			<< "-c chip" << endl
			<< "-g goenr" << endl
			<< "-r exprregenr" << endl
			<< "-d dnase" << endl;
		return 0;
	}
	Framework fw;
	if(fw.init(argc,argv)!=0)
	{
		return 0;
	}
	fw.genClusterAttribs();
	fw.genClusterRegulators();
	return 0;

}

