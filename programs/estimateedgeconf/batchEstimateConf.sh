if [ $# -ne 1 ] 
then
	echo "Usage batchEstimateConf outputdir"
	exit
fi
OUTPUTDIR=$1

#for DATA in lower upper lower_pool upper_pool
for DATA in ~/results/networkinference/gasch_causton/bootstrap/infnet_files
do
	NWFILE=$DATA.txt
	for THRESH in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 
	do
		RESDIR=$OUTPUTDIR/thresh$THRESH/$DATA
		echo "mkdir -p $RESDIR"
		mkdir -p $RESDIR
		echo "./estimateEdgeConf $NWFILE $THRESH $RESDIR/aggregate.txt filterededges"
		./estimateEdgeConf $NWFILE $THRESH $RESDIR/aggregate.txt filterededges
	done
done
