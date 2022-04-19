if [ $# -ne 1 ] 
then
	echo "Usage batchEstimateConf outputdir"
	exit
fi
OUTPUTDIR=$1
GOFILE=/home/local/MORGRIDGE/sroy/data/mouse/go/go_v2/geneontology_ensemblid_regnet.txt
GOFILE=/home/local/MORGRIDGE/sroy/data/mouse/mousegotermap_regnet.txt
HIGHCONFREGNET=/home/local/MORGRIDGE/sroy/results/networkinference/pgpm/rupa/r4_p-5_h0.6/model1/edgec_0.6_regnet.txt
HIGHCONFREGNET=/home/local/MORGRIDGE/sroy/results/networkinference/pgpm/influenza/analysis_mar14/withsignaling/edgec_1_regnet.txt
MSIGDB=/home/local/MORGRIDGE/sroy/data/msigdb/genenames_msigpath_gpl7202_regnet.txt
MSIGDBMOTIF=/home/local/MORGRIDGE/sroy/data/mouse/msigdb/mouse_motifs_regnet.txt
MOTIF=/home/local/MORGRIDGE/sroy/data/mouse/mouse_motifs_msigdb_jaspar_regnet.txt
CHIPNET=/home/local/MORGRIDGE/sroy/data/DCdata/gene_promoter_regnet.txt
DENSECLUSTERS=/home/local/MORGRIDGE/sroy/projects/mouse_attie_dave_coon/dense_clusters/biogrid_.7_density.gene_regnet.txt
#KONET=/home/local/MORGRIDGE/sroy/data/DCdata/konetout_1_tftgt_regnet.txt
KONET=/home/local/MORGRIDGE/sroy/data/DCdata/konetout_1_tftgt_regnet.txt
ITER=1

while [ $ITER -le 1 ]
do
	#RESDIR=$OUTPUTDIR/fold${ITER}
	#RESDIR=$OUTPUTDIR/randinit${ITER}/r4_p-5_h0.6/model1/fold0
	RESDIR=$OUTPUTDIR/randinit${ITER}/
	#OUTSUFF=$RESDIR/goproc
	#OUTSUFF=$RESDIR/tfenr.txt
	#GENELIST=$RESDIR/dc${ITER}_clusterassign.txt
	GENELIST=$RESDIR/mRNA.txt
	NWFILE=$RESDIR/mRNAsets.txt
	#echo "./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $RESDIR/goproc persg"
	#./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $RESDIR/mRNA_goproc persg
	GENELIST=$RESDIR/protein.txt
	NWFILE=$RESDIR/proteinsets.txt
	#./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $RESDIR/protein_goproc persg
	GENELIST=$RESDIR/phosphoprotein.txt
	NWFILE=$RESDIR/phosphoproteinsets.txt
	#./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $RESDIR/phosphoprotein_goproc persg
	GENELIST=$RESDIR/clusterassign_commonnames.txt
	NWFILE=$RESDIR/genesets_common.txt
	./enrichAnalyzer $NWFILE $GENELIST $DENSECLUSTERS 0.05 $RESDIR/densecluster_biogrid persg
	#./enrichAnalyzer $NWFILE $GENELIST $HIGHCONFREGNET 0.05 $RESDIR/module_reginfenr persg
	#./enrichAnalyzer $NWFILE $GENELIST $CHIPNET 0.05 $RESDIR/module_chipenr persg
	#./enrichAnalyzer $NWFILE $GENELIST $MSIGDB 0.05 $RESDIR/msigdb_enr persg
	#echo "./enrichAnalyzer $NWFILE $GENELIST $KONET 0.05 $RESDIR/module_koenr persg"
	#./enrichAnalyzer $NWFILE $GENELIST $KONET 0.05 $RESDIR/module_koenr persg
	#./enrichAnalyzer $NWFILE $GENELIST $MSIGDBMOTIF 0.05 $RESDIR/module_msigdb_motifenr persg
	#./enrichAnalyzer $NWFILE $GENELIST $MOTIF 0.05 $RESDIR/module_msigdb_jaspar_motifenr persg
	ITER=`expr $ITER + 1`
done

