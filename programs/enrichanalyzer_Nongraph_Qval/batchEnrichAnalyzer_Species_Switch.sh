if [ $# -ne 1 ] 
then
	echo "Usage batchEstimateConf outputdir"
	#exit
fi
#DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/carbonclusters_dup2c/allclade4
#DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/osrspecies/carbon
DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/carbon_mergedassign_wkp3/allclade4
#DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/schizodata_all4/soft_lf0.8_nlf0.8
DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/stressdata/ilan_heatstress_8spec_pp_postem/
GOFILE=/seq/compbio/sroy/networklearning/data/goslim/goslimprocgotermap.txt
GOFILE=/seq/compbio/sroy/networklearning/data/go/goproc
#GOFILE=/ahg/regev/users/sroy/compfuncgen/data/schizodata/pombegoloc_dash.txt
#SPECLIST=../scripts/speclist4.txt
SPECLIST=../scripts/speclist_heat_8spec.txt
#SPECLIST=../scripts/speclist_osr.txt
#SPECLIST=../scripts/speclist_schizo.txt
K=5
while [ $K -le 5 ]
do	
	for LEAF in  0.8
	do
	for RITER in 1 2 3 4 5
	do
	#for NONLEAF in 0.2 0.4 0.6 0.8
	#do
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
		NWFILE=$RESDIR/${SPECIES}_clusterswitch.txt
		OUTSUFF=$RESDIR/${SPECIES}_clusterswitch_goprocenr
		#OUTSUFF=$RESDIR/${SPECIES}_goslimprocenr
		GENELIST=$RESDIR/${SPECIES}_clusterassign.txt
		echo "bsub -o scratch/gmm_cc_osr_k$ITER.txt ./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $OUTSUFF persg"
		bsub -q hour -o scratch/gmm_allclade_k$ITER.txt ./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $OUTSUFF persg
	#done
	done
	done
	done
	K=`expr $K + 2`
done

