if [ $# -lt 1 ]
then
	echo "./batchReorder OUTDIR"
	exit
fi

REORDER=../reordermat/reorderMat

for CONS in 0.0 0.2 0.4 0.6 0.8
do
	for MOECNT in 10 15 20 25
	do
		ADJFNAME=$1/adj_cons${CONS}_k${MOECNT}.CDT
		ADJFNAME_MAT=$1/adj_cons${CONS}_k${MOECNT}_mat.txt
		echo "$REORDER $ADJFNAME $ADJFNAME_MAT"
		$REORDER $ADJFNAME $ADJFNAME_MAT
		RUN=`expr $RUN  + 1`
	done
done
