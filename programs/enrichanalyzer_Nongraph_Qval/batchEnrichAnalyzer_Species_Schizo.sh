if [ $# -ne 1 ] 
then
	echo "Usage batchEstimateConf outputdir"
	#exit
fi
DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/osrspecies/carbon
#DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/carbonclusters/allclade0
DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/schizodata_all4/merged_assign
#GOFILE=/seq/compbio/sroy/networklearning/data/goslim/goslimprocgotermap.txt
GOFILE=/ahg/regev/users/sroy/compfuncgen/data/schizodata/pombegoslimproc_dash.txt
#GOFILE=/ahg/regev/users/sroy/compfuncgen/data/schizodata/pombegoslimproc_inexact.txt
SPECLIST=../scripts/speclist0.txt
SPECLIST=../scripts/speclist_osr.txt
SPECLIST=../scripts/speclist_schizo.txt
K=8
while [ $K -le 8 ]
do
	RITER=3
	while [ $RITER -le 3 ]
	do
	RESDIR=$DIR/crossspecies_k${K}/randinit${RITER}
	for SPECIES in `cat $SPECLIST`
	do
		SUBSTR=Anc
		LEN=`expr match $SPECIES $SUBSTR`
		if [ $LEN -gt 0 ]
		then
			echo "Continuing at $SPECIES $LEN"
		#	continue
		fi
		NWFILE=$RESDIR/${SPECIES}_genesets.txt
		#OUTSUFF=$RESDIR/${SPECIES}_go_inexact_enr
		OUTSUFF=$RESDIR/${SPECIES}_goenr_strict
		GENELIST=$RESDIR/${SPECIES}_clusterassign.txt
		SCRATCHFILE=scratch/gmm_schizo_k${K}_r${RITER}.txt
		if [ -f $SCRATCHFILE ]
		then
			echo "rm $SCRATCHFILE"
			rm $SCRATCHFILE
		fi
		echo "bsub -o scratch/gmm_cc_osr_k$ITER.txt ./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $OUTSUFF persg"
		bsub -q hour -o $SCRATCHFILE ./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $OUTSUFF persg
	done
	RITER=`expr $RITER + 1`
	done
	K=`expr $K + 2`
done

