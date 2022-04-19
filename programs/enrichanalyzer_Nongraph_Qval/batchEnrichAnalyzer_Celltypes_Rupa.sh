if [ $# -ne 1 ] 
then
	echo "Usage batchEstimateConf outputdir"
	#exit
fi
DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/carbonclusters_dup2c/allclade4
#DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/osrspecies/carbon
DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_notcrosspecies/carbon_mergedassign_wkp5/allclade4
DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/carbon_mergedassign_wkp3/allclade4
DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/carbonclusters_dup2c/allclade4
#DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/schizodata_all4/soft_lf0.8_nlf0.8
#DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/carbon_ribospecies/
#DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/stressdata/ilan_heatstress_8spec_pp_postem/
#DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/stressdata/ilan_heatstress_8spec_calb42_cgla37_pp_postem/
#DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/stressdata/ilan_carbonstress_8spec_pp_postem/
#DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/PLEASE_ARBORETUM_ME/osrspec/
#DIR=/seq/compbio/sroy/networklearning/softclustering/clustering/ilandata_run/carbondata/
GOFILE=/seq/compbio/sroy/networklearning/data/goslim/goslimprocgotermap.txt
GOFILE=../../../data/human/humango/ensemblid_instance_go.txt
GOFILE=/work/sroy/kellis_lab/compbio/networklearning_old/data/go/goproc
#GOFILE=/seq/compbio/sroy/networklearning/data/go/riboprotogenes_regnet.txt
#GOFILE=/seq/compbio/sroy/networklearning/data/go/ribobiogenes_regnet.txt
#GOFILE=/ahg/regev/users/sroy/compfuncgen/data/schizodata/pombegoproc_dash.txt
SPECLIST=../scripts/speclist4.txt
DIR=../../../results/compfuncgen/cross_speciescluster/mammals/collapsed
DIR=../../../results/compfuncgen/cross_speciescluster/allstress/
DIR=../../../results/compfuncgen/cross_cellcluster/analysis_feb20_2014/
GOFILE=~/data/mouse/mousegotermap_regnet.txt
#SPECLIST=../scripts/speclist4_wkp.txt
#SPECLIST=../scripts/speclist_ribospecies.txt
#SPECLIST=../scripts/speclist_heat_8spec.txt
#SPECLIST=../scripts/speclist_osr.txt
#SPECLIST=../scripts/speclist_schizo.txt
K=20
while [ $K -le 20 ]
do	
	for RITER in  1
	do
	RESDIR=$DIR/k${K}/randinit${RITER}/
	for SPECIES in  ips pips mef A1 A2
	do
		SUBSTR=Anc
		LEN=`expr match $SPECIES $SUBSTR`
		if [ $LEN -gt 0 ]
		then
			echo "Continuing at $SPECIES $LEN"
		#	continue
		fi
		#NWFILE=$RESDIR/${SPECIES}_genesets_${SPECIES}.txt
		NWFILE=$RESDIR/${SPECIES}_genesets.txt
		OUTSUFF=$RESDIR/${SPECIES}_goprocenr
		GENELIST=$RESDIR/${SPECIES}_clusterassign.txt
		echo "./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $OUTSUFF persg"
		./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $OUTSUFF persg
	#done
	done
	done
	K=`expr $K + 5`
done

