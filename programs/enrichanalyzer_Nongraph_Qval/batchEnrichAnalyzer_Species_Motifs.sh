if [ $# -ne 1 ] 
then
	echo "Usage batchEstimateConf outputdir"
	#exit
fi
DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/osrspecies/carbon
DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/carbon_mergedassign_wkp5/allclade4
DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/carbonclusters_dup2c/allclade4
#DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/PLEASE_ARBORETUM_ME/osrspec/
#DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/stressdata/ilan_heatstress_8spec_calb42_cgla37_pp_postem/
#DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/schizodata_all4/soft_lf0.8_nlf0.8
#GOFILE=/seq/compbio/sroy/networklearning/data/goslim/goslimprocgotermap.txt
MOTIFDIR=../../../data/compfuncgen/jay_del_mutant/DGE/
#GOFILE=/ahg/regev/users/sroy/compfuncgen/data/schizodata/pombegoslimproc_dash.txt
#SPECLIST=../scripts/speclist4_wkp.txt
SPECLIST=../scripts/speclist4.txt
DIR=../../../results/compfuncgen/cross_speciescluster/osr_8spec/prefilter_zeromean
#SPECLIST=../scripts/speclist_osr.txt
#SPECLIST=../scripts/speclist_heat_8spec.txt
K=5
while [ $K -le 5 ]
do	
	for LEAF in  0.8
	do
	#for RITER in 1 2 3 4 5
	for RITER in  1 2 3 4 5
	do
	#for NONLEAF in 0.2 0.4 0.6 0.8
	#do
	
	RESDIR=$DIR/k${K}/randinit${RITER}/
	#RESDIR=$DIR/crossspecies_k${K}/randinit${RITER}/
	#for SPECIES in Scer Spom Cgla Scas Klac Calb
	for SPECIES in Scer Spom Cgla Scas Klac Calb
	do
		SUBSTR=Anc
		LEN=`expr match $SPECIES $SUBSTR`
		if [ $LEN -gt 0 ]
		then
			echo "Continuing at $SPECIES $LEN"
			continue
		fi
		LOWERSPEC=`echo $SPECIES | awk '{printf("%s\n",tolower($1))}'`
		MOTIFFILE=$MOTIFDIR/${LOWERSPEC}_regnet.txt
		OUTSUFF=$RESDIR/${SPECIES}_ko
		if [ $SPECIES == Scer ]
		then
			#GENELIST=$RESDIR/${SPECIES}_speciesspecnames_clusterassign.txt
			GENELIST=$RESDIR/${SPECIES}_speciesspecnames_clusterassign.txt
			NWFILE=$RESDIR/${SPECIES}_genesets.txt
		else
			#GENELIST=$RESDIR/${SPECIES}_speciesspecnames_clusterassign.txt
			GENELIST=$RESDIR/${SPECIES}_speciesspecnames_clusterassign.txt
			NWFILE=$RESDIR/${SPECIES}_genesets.txt
		fi
		echo "./enrichAnalyzer $NWFILE $GENELIST $MOTIFFILE 0.05 $OUTSUFF persg"
		./enrichAnalyzer $NWFILE $GENELIST $MOTIFFILE 0.05 $OUTSUFF persg
	#done
	done
	done
	done
	K=`expr $K + 2`
done

