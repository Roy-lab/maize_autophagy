SRC =    Framework.C  GeneNameMapper.C GeneExpManager.C optimalleaforder/HierarchicalCluster.C optimalleaforder/HierarchicalClusterNode.C optimalleaforder/OptimalLeafOrder.C optimalleaforder/Heap.C optimalleaforder/Distance.C

CC=g++
CFLAGS = -g
INCLFLAG = -I optimalleaforder

genClusterAttrib: $(SRC)
	$(CC) $(SRC)  $(CFLAGS) $(INCLFLAG) -o genClusterAttrib
clean: 
	rm -f genClusterAttrib
	rm -f *~
