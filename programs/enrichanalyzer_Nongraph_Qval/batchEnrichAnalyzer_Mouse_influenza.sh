if [ $# -ne 1 ] 
then
	echo "Usage batchEstimateConf outputdir"
	exit
fi
OUTPUTDIR=$1
GOFILE=~/data/mouse/mousegotermap_regnet.txt
HIGHCONFREGNET=/home/local/MORGRIDGE/sroy/results/networkinference/pgpm/influenza/analysis_mar14/withsignaling/edgec_1_regnet.txt
HIGHCONFREGNET=~/results/networkinference/pgpm/influenza/mouse_include_nonreg_influenza_d45/edgec_0.3_regnet.txt
MSIGDB=~/data/mouse/msigdb/c2.all.v4.0.symbols_filtered_mousenames_regnet.txt
MSIGDBMOTIF=~/data/mouse/msigdb/mouse_motifs_mappedtfs_v2_regnet.txt
MOTIF=~/data/mouse/mouse_motifs_msigdb_jaspar_regnet.txt
CHIPNET=~/data/DCdata/gene_promoter_regnet.txt
#KONET=/home/local/MORGRIDGE/sroy/data/DCdata/konetout_1_tftgt_regnet.txt
KONET=~/data/DCdata/konetout_1_tftgt_regnet.txt
SCREENNET=~/data/influenza/screen_gene_sets/mouse_screen_regnet.txt
ITER=0

#while [ $ITER -le 0 ]
#for ITER in c100 c160 
for ITER in h0.5sim h0.25sim h0.75sim
do
	#RESDIR=$OUTPUTDIR/fold${ITER}
	#RESDIR=$OUTPUTDIR/randinit${ITER}/r4_p-5_h0.6/model1/fold0
	#RESDIR=$OUTPUTDIR/part${ITER}/model1/fold0
	#RESDIR=$OUTPUTDIR/spectral_postproc/${ITER}
	RESDIR=$OUTPUTDIR/hierarchical_postproc/${ITER}
	#NWFILE=$RESDIR/genesets_c0.3.txt
	#NWFILE=$RESDIR/mouse_${ITER}_geneset.txt
	NWFILE=$RESDIR/geneset.txt
	#OUTSUFF=$RESDIR/goproc
	#OUTSUFF=$RESDIR/tfenr.txt
	#GENELIST=$RESDIR/dc${ITER}_clusterassign.txt
	#GENELIST=$RESDIR/modules_c0.3.txt
	GENELIST=$OUTPUTDIR/targets_c0.3.txt
	#HIGHCONFREGNET=$RESDIR/var_mb_regnet_prot.txt
	#./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $RESDIR/goproc persg
	./enrichAnalyzer $NWFILE $GENELIST $HIGHCONFREGNET 0.05 $RESDIR/module_reginfenr persg
	#./enrichAnalyzer $NWFILE $GENELIST $CHIPNET 0.05 $RESDIR/module_chipenr persg
	#./enrichAnalyzer $NWFILE $GENELIST $MSIGDB 0.05 $RESDIR/msigdb_enr persg
	echo "./enrichAnalyzer $NWFILE $GENELIST $KONET 0.05 $RESDIR/module_koenr persg"
	#./enrichAnalyzer $NWFILE $GENELIST $KONET 0.05 $RESDIR/module_koenr persg
	echo "./enrichAnalyzer $NWFILE $GENELIST $SCREENNET 0.05 $RESDIR/screenenr persg"
	#./enrichAnalyzer $NWFILE $GENELIST $SCREENNET 0.05 $RESDIR/screenenr persg
	#./enrichAnalyzer $NWFILE $GENELIST $MSIGDBMOTIF 0.05 $RESDIR/module_msigdb_motifenr persg
	#./enrichAnalyzer $NWFILE $GENELIST $MOTIF 0.05 $RESDIR/module_msigdb_jaspar_motifenr persg
	#ITER=`expr $ITER + 1`
done

