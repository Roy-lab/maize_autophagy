if [ $# -ne 1 ] 
then
	echo "Usage batchEstimateConf outputdir"
	exit
fi
OUTPUTDIR=$1
GOFILE=../../../data/human/humango/geneontology_cnames_regnet.txt
MSIGDB_CP=/mnt/ws/sysbio/roygroup/shared/data_new/human/genesets/c2.all.v4.0.symbols_filtered_regnet.txt
MSIGDB_MOTIF=/mnt/ws/sysbio/roygroup/shared/data_new/human/genesets/c3.all.v4.0.symbols_regnet.txt
MSIGDB_MOTIF=~/data/human/msigdb_motifs/msigdb_motifs_tfnames_mapped_regnet.txt
MOTIFNET=/mnt/ws/sysbio/roygroup/shared/data_new/human/motif/motif_promoter/transcript_based/human_motif_net_common_final_target_symbols_regnet.txt
SCREENNET=~/data/influenza/screen_gene_sets/human_screen_withflulist_regnet.txt
HIGHCONFREGNET=$OUTPUTDIR/edgec_0.3_regnet.txt
ITER=0

#while [ $ITER -le 0 ]
#for ITER in c100 c160
for ITER in h0.25sim h0.5sim h0.75sim
#for ITER in 0 1 2
do
	#RESDIR=$OUTPUTDIR/fold${ITER}
	#RESDIR=$OUTPUTDIR/part${ITER}/model1/fold0
	#RESDIR=$OUTPUTDIR/spectral_postproc/${ITER}
	RESDIR=$OUTPUTDIR/hierarchical_postproc/${ITER}
	#NWFILE=$RESDIR/genesets_c0.3.txt
	#NWFILE=$RESDIR/human_${ITER}_geneset.txt
	NWFILE=$RESDIR/geneset.txt
	GENELIST=$OUTPUTDIR/targets_c0.3.txt
	echo "./enrichAnalyzer $NWFILE $GENELIST $MOTIFNET 0.05 $RESDIR/motifenr persg"
	#./enrichAnalyzer $NWFILE $GENELIST $MOTIFNET 0.05 $RESDIR/motifenr persg
	echo "./enrichAnalyzer $NWFILE $GENELIST $MSIGDB_CP 0.05 $RESDIR/msigdbenr persg"
	./enrichAnalyzer $NWFILE $GENELIST $MSIGDB_CP 0.05 $RESDIR/msigdb_filtered_enr persg
	echo "./enrichAnalyzer $NWFILE $GENELIST $MSIGDB_MOTIF 0.05 $RESDIR/msigdb_motifenr persg"
	./enrichAnalyzer $NWFILE $GENELIST $MSIGDB_MOTIF 0.05 $RESDIR/msigdb_motifenr persg
	echo "./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $RESDIR/goenr persg"
	./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $RESDIR/goenr persg
	echo "./enrichAnalyzer $NWFILE $GENELIST $HIGHCONFREGNET 0.05 $RESDIR/modulereginfenr persg"
	./enrichAnalyzer $NWFILE $GENELIST $HIGHCONFREGNET 0.05 $RESDIR/modulereginfenr persg
	echo "./enrichAnalyzer $NWFILE $GENELIST $SCREEN 0.05 $RESDIR/screenenr persg"
	./enrichAnalyzer $NWFILE $GENELIST $SCREENNET 0.05 $RESDIR/screenenr persg
	#ITER=`expr $ITER + 1`
done

