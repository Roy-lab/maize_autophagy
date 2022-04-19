if [ $# -ne 1 ] 
then
	echo "Usage batchEstimateConf outputdir"
	exit
fi
OUTPUTDIR=$1

for DATA in lower upper lower_pool upper_pool
do
	NWFILE=$DATA.txt
	RESDIR=$OUTPUTDIR/$DATA
	THRESH=0
	echo "./estimateEdgeConf $NWFILE $THRESH $RESDIR/edgeconf alledges"
	./estimateEdgeConf $NWFILE $THRESH $RESDIR/edgeconf alledges
done
