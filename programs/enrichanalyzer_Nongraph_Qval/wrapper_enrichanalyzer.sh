#!/bin/sh
echo "Running here: $HOSTNAME"
tar -xvzf sharedlib.tgz
export LD_LIBRARY_PATH=sharedlib
echo "./enrichAnalyzer $1_genesets.txt $1_clusterassign.txt on_go_regnet.txt 0.05 ${1}_goenr persg"
./enrichAnalyzer $1_genesets.txt $1_clusterassign.txt on_go_regnet.txt 0.05 ${1}_goenr persg
