if [ $# -ne 1 ] 
then
	echo "Usage batchEstimateConf outputdir"
	exit
fi
OUTPUTDIR=$1
HIGHCONFREGNET=$OUTPUTDIR/edgec_0.4_regnet.txt
#MSIGDB=/home/local/MORGRIDGE/sroy/data/msigdb/genenames_msigpath_gpl7202_regnet.txt
MOTIFNET=../../../data/regnets/gordon_yetfasco_motifs_testmotif_union_regnet.txt
CHIPNET=../../../data/regnets/macisaac/macisaac_p0.001_cons0_regnet.txt
#CHIPNET=../../../data/regnets/macisaac/macisaac_p0.001_cons0_regnet_met4_32.txt
#GOFILE=/mnt/ws/sysbio/roygroup/shared/data_new/yeast/annotations/GO/goproc
#KONET=/home/local/MORGRIDGE/sroy/data/DCdata/konetout_1_tftgt_regnet.txt
#KONET=/home/local/MORGRIDGE/sroy/data/DCdata/konetout_1_tftgt_regnet.txt
ITER=0

while [ $ITER -le 4 ]
do
	#RESDIR=$OUTPUTDIR/fold${ITER}
	#RESDIR=$OUTPUTDIR/randinit${ITER}/r4_p-5_h0.6/model1/fold0
	RESDIR=$OUTPUTDIR/
	RESDIR=$OUTPUTDIR//part${ITER}/model1/fold${ITER}
	#NWFILE=$RESDIR/module_regulator${ITER}_geneset.txt
	NWFILE=$RESDIR/genesets_highconf_c0.4.txt
	#GENELIST=$RESDIR/dc${ITER}_clusterassign.txt
	GENELIST=$OUTPUTDIR/highconf_targets.txt
	#GENELIST=$OUTPUTDIR/r4_p-5_h0.6/module_regulator${ITER}_geneset.txt
	#GENELIST=$RESDIR/segal_yeast_fold${ITER}_clusterassign.txt
	#./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $OUTSUFF persg
	#./enrichAnalyzer $NWFILE $GENELIST $HIGHCONFREGNET 0.05 $RESDIR/module_reginfenr persg
	#./enrichAnalyzer $NWFILE $GENELIST $CHIPNET 0.05 $RESDIR/module_chipenr persg
	#echo "./enrichAnalyzer $NWFILE $GENELIST $MOTIFNET 0.05 $OUTSUFF persg"
	#./enrichAnalyzer $NWFILE $GENELIST $MOTIFNET 0.05 $OUTSUFF persg
	#OUTSUFF=$RESDIR/chipenr_met4_met32
	OUTSUFF=$RESDIR/chipenr_module_restrict${ITER}
	echo "./enrichAnalyzer $NWFILE $GENELIST $CHIPNET 0.05 $OUTSUFF persg"
	#./enrichAnalyzer $NWFILE $GENELIST $CHIPNET 0.05 $OUTSUFF persg
	OUTSUFF=$RESDIR/motifenr
#	./enrichAnalyzer $NWFILE $GENELIST $MOTIFNET 0.05 $OUTSUFF persg
	OUTSUFF=$RESDIR/goproc
#	./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $OUTSUFF persg
	#OUTSUFF=$RESDIR/module_reginfenr${ITER}
	OUTSUFF=$RESDIR/module_reginfenr
	./enrichAnalyzer $NWFILE $GENELIST $HIGHCONFREGNET 0.05 $OUTSUFF persg
	#OUTSUFF=$RESDIR/module_allreginfenr
	#REGNET=$RESDIR/var_mb_pw_k299.txt
	#awk '{printf("%s\t%s\n",$2,$1)}' $REGNET > regnet.txt
	#./enrichAnalyzer $NWFILE $GENELIST regnet.txt 0.05 $OUTSUFF persg
	#rm regnet.txt
	#./enrichAnalyzer $NWFILE $GENELIST $MSIGDB 0.05 $RESDIR/msigdb_enr persg
	#echo "./enrichAnalyzer $NWFILE $GENELIST $KONET 0.05 $RESDIR/module_koenr persg"
	ITER=`expr $ITER + 1`
done

