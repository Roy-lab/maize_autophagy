if [ $# -ne 1 ] 
then
	echo "Usage batchEstimateConf outputdir"
	#exit
fi
DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/osrspecies/carbon
DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/carbon_mergedassign_wkp3/allclade4
DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/carbonclusters_dup2c/allclade4
#DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/schizodata_all4/soft_lf0.8_nlf0.8
GOFILE=/seq/compbio/sroy/networklearning/data/go/goproc
#MOTIFFILE=/ahg/regev/users/sroy/compfuncgen/data/schizodata/pombegoslimproc_dash.txt
#SPECLIST=../scripts/speclist4_wkp.txt
SPECLIST=../scripts/speclist4.txt
#SPECLIST=../scripts/speclist_osr.txt
#SPECLIST=../scripts/speclist_schizo.txt
K=5
while [ $K -le 5 ]
do	
	for LEAF in  0.8
	do
	for RITER in 4
	do
	RESDIR=$DIR/crossspecies_k${K}/randinit${RITER}/lf${LEAF}_nlf${LEAF}
	for SPECIES in `cat $SPECLIST`
	do
		SUBSTR=Anc
		LEN=`expr match $SPECIES $SUBSTR`
		if [ $LEN -gt 0 ]
		then
			echo "Continuing at $SPECIES $LEN"
		#	continue
		fi
		OUTSUFF=$RESDIR/${SPECIES}_dup_goenr
		if [ $SPECIES == Scer ]
		then
			GENELIST=$RESDIR/duplistswitchbg_${SPECIES}.txt
			NWFILE=$RESDIR/duplist_${SPECIES}.txt
		else
			GENELIST=$RESDIR/duplistswitchbg_${SPECIES}.txt
			NWFILE=$RESDIR/duplist_${SPECIES}.txt
		fi
		echo "bsub -q hour -o scratch/gmm_allclade_${SPECIES}_k${ITER}.txt ./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $OUTSUFF persg"
		bsub -q hour -o scratch/gmm_allclade_${SPECIES}_k$ITER.txt ./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $OUTSUFF persg
	done
	done
	done
	K=`expr $K + 2`
done

