if [ $# -ne 1 ] 
then
	echo "Usage batchEstimateConf outputdir"
#	exit
fi
OUTPUTDIR=/seq/compbio/sroy/btomics/jonathan/analysis/gmmclusters/unnormalized/
GOFILE=/seq/compbio/sroy/btomics/jonathan/data/cog/cog_goslim.txt

for K in 25 35 45
do
ITER=1
while [ $ITER -le 5 ]
do
	RESDIR=$OUTPUTDIR/k${K}/randinit${ITER}
	NWFILE=$RESDIR/genesets.txt
	OUTSUFF=$RESDIR/cogenr.txt
	GENELIST=$RESDIR/clusterassign_zvalues.txt
	SCRATCH=scratch/bt_${ITER}.txt
	if [ -f $SCRATCH ]
	then 
		rm $SCRATCH
	fi
	echo "bsub -o $SCRATCH ./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $OUTSUFF persg"
	bsub -o $SCRATCH ./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $OUTSUFF persg
	ITER=`expr $ITER + 1`
done
done
