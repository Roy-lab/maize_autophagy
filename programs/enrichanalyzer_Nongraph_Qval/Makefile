LFLAG = -lgsl -lgslcblas 
SRC =  Clique.C  Distance.C  Framework.C  GeneNameMapper.C  GOMgr.C  GOTerm.C	Graph.C  HyperGeomPval.C  Randomizer.C	SubGraph.C  SubGraphMgr.C  Vertex.C 
LIBPATH = /ahg/regev/users/sroy/thirdparty/gsl_install/lib
INCLPATH1 = /ahg/regev/users/sroy/thirdparty/gsl_install/include

CC=g++
CFLAGS = -g

enrichAnalyzer: $(SRC)
	$(CC) $(SRC) -I $(INCLPATH1)  -L $(LIBPATH) $(LFLAG) $(CFLAGS) -o enrichAnalyzer
	cp enrichAnalyzer ../../
clean: 
	rm enrichAnalyzer
