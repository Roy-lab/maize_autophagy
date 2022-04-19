if [ $# -ne 1 ] 
then
	echo "Usage batchEstimateConf outputdir"
	exit
fi
OUTPUTDIR=$1
GOFILE=/seq/compbio/sroy/networklearning/data/goslim/goslimprocgotermap.txt
TFFILE=/seq/compbio/sroy/networklearning/data/tfnet/regnet_enr.txt
TFFILE=/seq/compbio/sroy/networklearning/data/combinedmotifs_ilanjay/combinedmotifs_regnet_enr.txt
#RANDFILE=/seq/compbio/sroy/networklearning/data/randgo/randgoproc_s3-44_n1000_k20.txt
#GENELIST=/ahg/regev/users/sroy/compfuncgen/data/gaschstress/subdatasets/heatshock/heatshock.genelist
#GENELIST=/ahg/regev/users/sroy/compfuncgen/data/stressdata/jaydata/species_specific/alldata/OSstress/scer/scer_high2_std2.genelist
GENELIST=/seq/compbio/sroy/networklearning/fgraph/data/cellcycle/cellcycle_30exp_high2_std2.genelist
ITER=1

for BETA1 in 0 0.5 1
do
	for BETA2 in 0  1.6  3.2
	do
	RESDIR=$OUTPUTDIR/k1/beta${BETA1}_${BETA2}/model1/fold0/
	NWFILE=$RESDIR/genesets.txt
	OUTSUFF=$RESDIR/combmotif
	SCRATCHFILE=scratch/cc_combmotis_$ITER.txt
	if [ -f $SCRATCHFILE ]
	then
		rm $SCRATCHFILE
	fi
	echo "bsub -o $SCRATCHFILE ./enrichAnalyzer $NWFILE $GENELIST $TFFILE 0.05 1e-5 2 0.05 $OUTSUFF"
	bsub -o $SCRATCHFILE ./enrichAnalyzer $NWFILE $GENELIST $TFFILE 0.05 1e-5 2 0.05 $OUTSUFF
	ITER=`expr $ITER + 1`
	done
done

