#for N in 200 300 400 500
for N in 1000
do
INPUTMODULES=../../../data/genenetweaver/p0.4/net$N/net${N}_geneset.txt
GENIEDIR=../../../results/networkinference/genie3/gnw/p0.4/net${N}/
for TOP in 0.2 0.3 0.4
do
GENIENET=$GENIEDIR/net${N}_Ntop${TOP}_regnet.txt
awk '{printf("%s\t%s\n",$2,$1)}' $GENIEDIR/net${N}_Ntop${TOP}.txt > $GENIENET
GENELIST=../../../data/genenetweaver/p0.4/net$N/net${N}_nw_clusters.txt
GENIEOUTSUFF=$GENIEDIR/truemodule_top${TOP}_reginfenr
echo "./enrichAnalyzer $INPUTMODULES $GENELIST $GENIENET 0.05 $GENIEOUTSUFF persg"
./enrichAnalyzer $INPUTMODULES $GENELIST $GENIENET  0.05 $GENIEOUTSUFF persg
done
done
