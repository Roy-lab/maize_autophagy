#for N in 100 200 300 400 500
for N in 1000
do
INPUTMODULES=../../../data/genenetweaver/p0.4/net$N/net${N}_geneset.txt
#PGPMDIR=../../../results/networkinference/pgpm/gnw/p0.4/net${N}/randpartitions/r8_p-1_h0.6/
PGPMDIR=../../../results/networkinference/pgpm/gnw/p0.4/net${N}/fivefold/r8_p-1_h0.6/
#PERGENEDIR=../../../results/networkinference/pergene/gnw/p0.4/net${N}/
PERGENEDIR=../../../results/networkinference/pergene/gnw/p0.4/net${N}/fivefold/p-1
MODNETDIR=/mnt/ws/sysbio/roygroup/shared/projects/softModuleNetworks/genenetweaver/p0.4/net${N}/results/
MODNETNET=net${N}_graphTau3_regnet.txt
awk '{printf("%s\t%s\n",$2,$1)}' $MODNETDIR/net${N}_graphTau3.txt > $MODNETNET
PGPMNET=$PGPMDIR/edgec_0.6_regnet.txt
awk '{printf("%s\t%s\n",$2,$1)}' $PGPMDIR/edge_c0.6.txt > $PGPMNET
PERGENENET=$PERGENEDIR/edgec_0.6_regnet.txt
awk '{printf("%s\t%s\n",$2,$1)}' $PERGENEDIR/edgec_0.6.txt > $PERGENENET


#GENIEDIR=../../../results/networkinference/genie3/gnw/p0.4/net${N}/
#GENIENET=$GENIEDIR/net${N}_Nfrompgpm_regnet.txt
#awk '{printf("%s\t%s\n",$2,$1)}' $GENIEDIR/net${N}_Nfrompgpm_r8.txt > $GENIENET
#GENIENET2=$GENIEDIR/net${N}_Nfrompergene_regnet.txt
#awk '{printf("%s\t%s\n",$2,$1)}' $GENIEDIR/net${N}_Nfrompergene.txt > $GENIENET2
GENELIST=../../../data/genenetweaver/p0.4/net$N/net${N}_nw_clusters.txt
PGPMOUTSUFF=$PGPMDIR/truemodule_reginfenr
PERGENEOUTSUFF=$PERGENEDIR/truemodule_reginfenr
MODNETOUTSUFF=$MODNETDIR/truemodule_reginfenr
#GENIEOUTSUFF=$GENIEDIR/truemodule_pgpm_reginfenr
#GENIEOUTSUFF2=$GENIEDIR/truemodule_pergene_reginfenr
echo "./enrichAnalyzer $INPUTMODULES $GENELIST $PGPMNET  0.05 $PGPMOUTSUFF persg"
./enrichAnalyzer $INPUTMODULES $GENELIST $PGPMNET  0.05 $PGPMOUTSUFF persg
echo "./enrichAnalyzer $INPUTMODULES $GENELIST $PERGENENET  0.05 $PERGENEOUTSUFF persg"
#./enrichAnalyzer $INPUTMODULES $GENELIST $PERGENENET  0.05 $PERGENEOUTSUFF persg
echo "./enrichAnalyzer $INPUTMODULES $GENELIST $MODNETNET  0.05 $MODNETOUTSUFF persg"
#./enrichAnalyzer $INPUTMODULES $GENELIST $MODNETNET  0.05 $MODNETOUTSUFF persg
#echo "./enrichAnalyzer $INPUTMODULES $GENELIST $GENIENET 0.05 $GENIEOUTSUFF persg"
#./enrichAnalyzer $INPUTMODULES $GENELIST $GENIENET 0.05 $GENIEOUTSUFF persg
#echo "./enrichAnalyzer $INPUTMODULES $GENELIST $GENIENET2 0.05 $GENIEOUTSUFF2 persg"
#./enrichAnalyzer $INPUTMODULES $GENELIST $GENIENET2 0.05 $GENIEOUTSUFF2 persg
done
