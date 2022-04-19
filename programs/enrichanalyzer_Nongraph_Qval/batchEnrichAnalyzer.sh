if [ $# -ne 1 ] 
then
	echo "Usage batchEstimateConf outputdir"
	exit
fi
OUTPUTDIR=$1
GOFILE=/seq/compbio/sroy/networklearning/data/goslim/goslimprocgotermap.txt
GOFILE=/seq/compbio/sroy/networklearning/data/go/goproc
TFFILE=/seq/compbio/sroy/networklearning/data/tfnet/regnet_enr.txt
TFFILE=/seq/compbio/sroy/networklearning/data/combinedmotifs_ilanjay/combinedmotifs_regnet_enr.txt
#RANDFILE=/seq/compbio/sroy/networklearning/data/randgo/randgoproc_s3-44_n1000_k20.txt
ITER=1

while [ $ITER -le 3 ]
do
	RESDIR=$OUTPUTDIR/randinit${ITER}
	NWFILE=$RESDIR/genesets.txt
	OUTSUFF=$RESDIR/goprocenr.txt
	#OUTSUFF=$RESDIR/tfenr.txt
	GENELIST=$RESDIR/clusterassign.txt
	SCRATCH=scratch/gs_hs_lambda6${ITER}.txt
	if [ -f $SCRATCH ]
	then 
		rm $SCRATCH
	fi
	echo "bsub -o scratch/gs_hs_lambda6$ITER.txt ./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $OUTSUFF persg"
	bsub -o scratch/gs_hs_$ITER.txt ./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $OUTSUFF persg
	ITER=`expr $ITER + 1`
done

