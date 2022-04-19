if [ $# -ne 1 ] 
then
	echo "Usage batchEstimateConf outputdir"
	exit
fi
OUTPUTDIR=$1
GOFILE=../../../data/plant/go/plant_go_regnet.txt
MOTIFNET=../../../data/plant/allplantmotif/known_motif_regnet.txt
HIGHCONFREGNET=$OUTPUTDIR/edgec_0.3_regnet.txt
ITER=0

while [ $ITER -le 4 ]
do
	#RESDIR=$OUTPUTDIR/fold${ITER}
	RESDIR=$OUTPUTDIR/part${ITER}/model1/fold0
	NWFILE=$RESDIR/genesets_c0.3.txt
	GENELIST=$OUTPUTDIR/targets_c0.3.txt
	echo "./enrichAnalyzer $NWFILE $GENELIST $MOTIFNET 0.05 $RESDIR/motifenr persg"
	./enrichAnalyzer $NWFILE $GENELIST $MOTIFNET 0.05 $RESDIR/motifenr persg
	echo "./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $RESDIR/goenr persg"
	./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $RESDIR/goenr persg
	echo "./enrichAnalyzer $NWFILE $GENELIST $HIGHCONFREGNET 0.05 $RESDIR/modulereginfenr persg"
	./enrichAnalyzer $NWFILE $GENELIST $HIGHCONFREGNET 0.05 $RESDIR/modulereginfenr persg
	ITER=`expr $ITER + 1`
done

