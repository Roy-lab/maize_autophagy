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
#MOTIFDIR=/home/local/MORGRIDGE/sroy/data/cichilddata/motifs/speciesspecific_regnets/
MOTIFDIR=~/data/cichilddata/motifs/refined_motifs/
SPECLIST=../scripts/speclist4.txt
DIR=../../../results/compfuncgen/cross_speciescluster/mammals/collapsed
DIR=../../../results/compfuncgen/cross_speciescluster/allstress/
DIR=../../../results/compfuncgen/cross_speciescluster/cichlids/
#SPECLIST=../scripts/speclist4_wkp.txt
#SPECLIST=../scripts/speclist_ribospecies.txt
#SPECLIST=../scripts/speclist_heat_8spec.txt
#SPECLIST=../scripts/speclist_osr.txt
#SPECLIST=../scripts/speclist_schizo.txt
K=15
while [ $K -le 15 ]
do	
	for RITER in  5 
	do
	#for NONLEAF in 0.2 0.4 0.6 0.8
	#do
	RESDIR=$DIR/k${K}/randinit${RITER}/
	OUTDIR=$RESDIR/refinedmotifs/
	mkdir $OUTDIR
	#RESDIR=$DIR/clusteroutput_k${K}/randinit${RITER}/iter2
	#RESDIR=$DIR/crossspecies_k${K}/randinit${RITER}/
	#for SPECIES in `cat $SPECLIST`
	#for SPECIES in  Human Chimpanzee Gorilla Macaque Mouse Opossum Platypus Chicken
	#for SPECIES in  Scer Cgla Scas Klac Calb Anc4 Anc5 Anc9 Anc11
	for SPECIES in  Ga Mz Nb On Pn Ab Anc1 Anc2 Anc3 Anc4 Anc5
	do
		SUBSTR=Anc
		LEN=`expr match $SPECIES $SUBSTR`
		if [ $LEN -gt 0 ]
		then
			echo "Continuing at $SPECIES $LEN"
			continue
		fi
		#NWFILE=$RESDIR/${SPECIES}_genesets_${SPECIES}.txt
		NWFILE=$RESDIR/${SPECIES}_speciesspecnames_genesets.txt
		OUTSUFF=$OUTDIR/${SPECIES}_motifenr
		GENELIST=$RESDIR/${SPECIES}_speciesspecnames_clusterassign.txt
		MOTIFFILE=$MOTIFDIR/${SPECIES}_motifnames_regnet.txt
		echo "./enrichAnalyzer $NWFILE $GENELIST $MOTIFFILE 0.05 $OUTSUFF persg"
		./enrichAnalyzer $NWFILE $GENELIST $MOTIFFILE 0.05 $OUTSUFF persg
	#done
	done
	done
	K=`expr $K + 2`
done

