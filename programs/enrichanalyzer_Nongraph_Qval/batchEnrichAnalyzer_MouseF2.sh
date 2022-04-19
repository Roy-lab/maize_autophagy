if [ $# -ne 1 ] 
then
	echo "Usage batchEstimateConf outputdir"
	exit
fi
OUTPUTDIR=$1
GOFILE=/home/local/MORGRIDGE/sroy/data/mouse/mousegotermap_regnet.txt
HIGHCONFREGNET=/home/local/MORGRIDGE/sroy/results/networkinference/pgpm/rupa/r4_p-5_h0.6/model1/edgec_0.6_regnet.txt
MSIGDB=/home/local/MORGRIDGE/sroy/data/msigdb/genenames_msigpath_gpl7202_regnet.txt
MSIGDBMOTIF=/home/local/MORGRIDGE/sroy/data/mouse/msigdb/mouse_motifs_regnet.txt
MOTIF=/home/local/MORGRIDGE/sroy/data/mouse/mouse_motifs_msigdb_jaspar_regnet.txt
CHIPNET=/home/local/MORGRIDGE/sroy/data/DCdata/gene_promoter_regnet.txt
#KONET=/home/local/MORGRIDGE/sroy/data/DCdata/konetout_1_tftgt_regnet.txt
KONET=/home/local/MORGRIDGE/sroy/data/DCdata/konetout_1_tftgt_regnet.txt
HIGHCONFREGNET=/home/local/MORGRIDGE/sroy/results/networkinference/pgpm/mouse_f2/liver/edgec_0.6_regnet.txt
ITER=0

while [ $ITER -le 5 ]
do
	#RESDIR=$OUTPUTDIR/fold${ITER}
	RESDIR=$OUTPUTDIR/randpartitions/part${ITER}/model1/fold0
	NWFILE=$RESDIR/genesets_highconf_c0.6.txt
	GENELIST=$RESDIR/modules_highconf_c0.6.txt
	./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $RESDIR/goproc persg
	./enrichAnalyzer $NWFILE $GENELIST $HIGHCONFREGNET 0.05 $RESDIR/module_reginfenr persg
	./enrichAnalyzer $NWFILE $GENELIST $MSIGDB 0.05 $RESDIR/msigdb_enr persg
	#echo "./enrichAnalyzer $NWFILE $GENELIST $KONET 0.05 $RESDIR/module_koenr persg"
	#./enrichAnalyzer $NWFILE $GENELIST $KONET 0.05 $RESDIR/module_koenr persg
	#./enrichAnalyzer $NWFILE $GENELIST $MSIGDBMOTIF 0.05 $RESDIR/module_msigdb_motifenr persg
	./enrichAnalyzer $NWFILE $GENELIST $MOTIF 0.05 $RESDIR/module_msigdb_jaspar_motifenr persg
	ITER=`expr $ITER + 1`
done

