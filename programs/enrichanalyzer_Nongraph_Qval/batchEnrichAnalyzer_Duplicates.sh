
if [ $# -ne 1 ] 
then
	echo "Usage batchEstimateConf outputdir"
	#exit
fi
DIR=/seq/compbio/sroy/networklearning/graphvalidation/estGeneralDupEnr/randinit4
GOFILE=/seq/compbio/sroy/networklearning/data/go/goproc
GOFILE=/seq/compbio/sroy/networklearning/data/motifs_t0.6/Scer_regnet.txt
for DUPSPECIES in Anc4 Anc5 Anc9 Anc10 Anc11 Ylip Anc13  all speciespec
#for DUPSPECIES in speciesspec
do
	NWFILE=$DIR/dupenr_${DUPSPECIES}_genes_switched_dup.txt
	#BGFILE=$DIR/dupenr_${DUPSPECIES}_genes_switched.txt
	BGFILE=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/carbonclusters_dup2c/allclade4/crossspecies_k5/randinit4/lf0.8_nlf0.8/Scer_clusterassign.txt
	OUTSUFF=$DIR/dupenr_${DUPSPECIES}_switched_motifenr_scer
	echo "./enrichAnalyzer $NWFILE $BGFILE $GOFILE 0.05 $OUTSUFF persg"
	./enrichAnalyzer $NWFILE $BGFILE $GOFILE 0.05 $OUTSUFF persg
	BGFILE=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/carbonclusters_dup2c/allclade4/crossspecies_k5/randinit4/lf0.8_nlf0.8/Scer_clusterassign.txt
	OUTSUFF=$DIR/dupenr_${DUPSPECIES}_notswitched_motifenr_scer
	NWFILE=$DIR/dupenr_${DUPSPECIES}_genes_notswitched_dup.txt
	./enrichAnalyzer $NWFILE $BGFILE $GOFILE 0.05 $OUTSUFF persg
done	
