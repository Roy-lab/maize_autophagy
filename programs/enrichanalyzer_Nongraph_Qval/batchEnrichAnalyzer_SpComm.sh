if [ $# -ne 1 ] 
then
	echo "Usage batchEstimateConf outputdir"
	exit
fi
OUTPUTDIR=$1
GOFILE=/seq/compbio/sroy/networklearning/data/go/goproc
#TFFILE=/seq/compbio/sroy/networklearning/data/tfnet/regnet_enr.txt
RANDFILE=/seq/compbio/sroy/networklearning/data/randgo/randgoproc_s3-44_n1000_k20.txt
ITER=1
for DATA in lucommon lupoolcommon lower_specific upper_specific lower_pool_specific upper_pool_specific
do
	for THRESH in 0.4 
	do
		RESDIR=$OUTPUTDIR/thresh$THRESH/$DATA
		for FSUFF in flat 
		do
			NWFILE=$RESDIR/targets_${FSUFF}.txt
			OUTSUFF=$RESDIR/targets_${FSUFF}_goproc
			FDRSUFF=$RESDIR/targets_${FSUFF}_fdr_goproc
			echo "bsub -o scratch/temp$ITER.txt ./enrichAnalyzer $NWFILE $OUTSUFF.txt 0.05 1e-5 2 $RANDFILE $GOFILE 0.05 $FDRSUFF"
			bsub -o scratch/temp$ITER.txt ./enrichAnalyzer $NWFILE $OUTSUFF.txt 0.05 1e-5 2 $RANDFILE $GOFILE 0.05 $FDRSUFF
			ITER=`expr $ITER + 1`
		done
	done
done

