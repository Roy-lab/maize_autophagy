#for N in 100 200 300 400 500
for N in 1000
do
INDIR=../../../data/genenetweaver/p0.4/net$N/
INPUTMODULES=$INDIR/net${N}_geneset.txt
GENELIST=../../../data/genenetweaver/p0.4/net$N/net${N}_nw_clusters.txt
INNET=$INDIR/net${N}_nw_nw.tsv
awk '{printf("%s\t%s\n",$2,$1)}' $INNET > $INDIR/net${N}_regnet.txt
OUTSUFF=$INDIR/processed/net${N}_regenr
echo "./enrichAnalyzer $INPUTMODULES $GENELIST $INDIR/net${N}_regnet.txt 0.05 $OUTSUFF persg"
./enrichAnalyzer $INPUTMODULES $GENELIST $INDIR/net${N}_regnet.txt  0.05 $OUTSUFF persg
awk '{printf("%s\t%s\n",$1,$2)}' ${OUTSUFF}_details.txt > $INDIR/processed/net${N}_regenr_nw.txt
done
