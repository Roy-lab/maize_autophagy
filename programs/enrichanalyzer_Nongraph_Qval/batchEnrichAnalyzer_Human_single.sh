if [ $# -ne 1 ] 
then
	echo "Usage batchEstimateConf outputdir"
	exit
fi
OUTPUTDIR=$1
GOFILE=../../../data/human/humango/geneontology_cnames_regnet.txt
#HIGHCONFREGNET=/home/local/MORGRIDGE/sroy/results/networkinference/pgpm/humanES/logratio/edgec_1_regnet.txt
HIGHCONFREGNET=/home/local/MORGRIDGE/sroy/results/networkinference/pgpm/humanES/zeromean/edgec_1_regnet.txt
#MSIGDB=/home/local/MORGRIDGE/sroy/data/msigdb/genenames_msigpath_gpl7202_regnet.txt
CHIPNET=/home/local/MORGRIDGE/sroy/data/human/tfnet/human_chipseq_regnet.txt
CHIPNET=/home/local/MORGRIDGE/sroy/data/human/tfnet/humanes_chip/H1_chip_regnet.txt
#DNASENET1=/home/local/MORGRIDGE/sroy/data/human/tfnet/dnase1/dnase1_thurman/Process1_HESC_regnet.txt
DNASENET2=/home/local/MORGRIDGE/sroy/data/human/tfnet/dnase1/dnase1_thurman/Process2_HESC_regnet.txt
DNASENET4=/home/local/MORGRIDGE/sroy/data/human/tfnet/dnase1/dnase1_thurman/Process4_HESC_regnet.txt
#KONET=/home/local/MORGRIDGE/sroy/data/DCdata/konetout_1_tftgt_regnet.txt
#KONET=/home/local/MORGRIDGE/sroy/data/DCdata/konetout_1_tftgt_regnet.txt
ITER=1

while [ $ITER -le 5 ]
do
	#RESDIR=$OUTPUTDIR/fold${ITER}
	#RESDIR=$OUTPUTDIR/randinit${ITER}/r4_p-5_h0.6/model1/fold0
	RESDIR=$OUTPUTDIR/
	NWFILE=$RESDIR/geneset.txt
	#OUTSUFF=$RESDIR/tfenr.txt
	GENELIST=$RESDIR/clusterassign.txt
	#GENELIST=$OUTPUTDIR/highconf_targets.txt
	OUTSUFF=$RESDIR/goproc
	./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $OUTSUFF persg
	#OUTSUFF=$RESDIR/module_reginfenr
	#./enrichAnalyzer $NWFILE $GENELIST $HIGHCONFREGNET 0.05 $OUTSUFF persg
	echo "./enrichAnalyzer $NWFILE $GENELIST $CHIPNET 0.05 $RESDIR/module_chipenr persg"
	./enrichAnalyzer $NWFILE $GENELIST $CHIPNET 0.05 $RESDIR/module_chipenr persg
	#./enrichAnalyzer $NWFILE $GENELIST $DNASENET1 0.05 $RESDIR/module_dnaseenr_proc1 persg
	#./enrichAnalyzer $NWFILE $GENELIST $DNASENET2 0.05 $RESDIR/module_dnaseenr_proc2 persg
	./enrichAnalyzer $NWFILE $GENELIST $DNASENET4 0.05 $RESDIR/module_dnaseenr_proc4 persg
	#./enrichAnalyzer $NWFILE $GENELIST $MSIGDB 0.05 $RESDIR/msigdb_enr persg
	#echo "./enrichAnalyzer $NWFILE $GENELIST $KONET 0.05 $RESDIR/module_koenr persg"
	ITER=`expr $ITER + 1`
done

