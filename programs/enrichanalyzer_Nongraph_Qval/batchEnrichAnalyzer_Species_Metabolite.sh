
if [ $# -ne 1 ] 
then
	echo "Usage batchEstimateConf outputdir"
	#exit
fi
#DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/metmad/data_10_22/tweakedsceronly_postprocessed_tweakedinit
DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/metmad/data_10_22/tweaked_pp_postem_lca
GOFILE=/ahg/regev/users/sroy/compfuncgen/data/metabolitemadness/metaboliteClasses_slimterm_Scer.txt
#GOFILE=/ahg/regev/users/sroy/compfuncgen/data/metabolitemadness/metpathways_slimterm_Scer.txt
SPECLIST=../scripts/speclist_metmad.txt
K=3
while [ $K -le 6 ]
do	
	for RITER in 1 2 3 4 5
	do
	RESDIR=$DIR/crossspecies_k${K}/randinit${RITER}/
	for SPECIES in `cat $SPECLIST`
	do
		SUBSTR=Anc
		LEN=`expr match $SPECIES $SUBSTR`
		if [ $LEN -gt 0 ]
		then
			echo "Continuing at $SPECIES $LEN"
			#continue
		fi
		NWFILE=$RESDIR/${SPECIES}_genesets.txt
		OUTSUFF=$RESDIR/${SPECIES}_metclassenr
		#OUTSUFF=$RESDIR/${SPECIES}_metpathenr
		GENELIST=$RESDIR/${SPECIES}_clusterassign.txt
		echo "bsub -o scratch/gmm_cc_osr_k$ITER.txt ./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $OUTSUFF persg"
		SCRATCHFILE=scratch/gmm_allclade_k$ITER.txt
		if [ -f $SCRATCHFILE ]
		then
			rm $SCRATCHFILE
		fi
		bsub -q hour -o $SCRATCHFILE ./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $OUTSUFF persg
		#./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $OUTSUFF persg
	done
	done
	K=`expr $K + 1`
done

