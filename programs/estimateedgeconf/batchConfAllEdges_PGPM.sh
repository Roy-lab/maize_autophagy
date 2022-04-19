if [ $# -ne 1 ] 
then
	echo "Usage batchEstimateConf outputdir"
	exit
fi
OUTPUTDIR=$1

#for H in 0.2 0.3 0.4 0.5 0.6 0.7 0.8
#do
#for P in -1 -2 -3 -4 -5 -6 -7 -8 -9 -10 -11 -12 -13 -14 -15 -20
#for P in -1  -5  -10 -15 -20
for P in -0.2 -0.4 -0.6 -0.8 -1
do
#for R in 2 4 8 10 12 14
#do
	#DIRNAME=r${R}_p${P}_h${H}
	DIRNAME=p${P}
	RESDIR=$OUTPUTDIR/$DIRNAME
	echo "find $RESDIR -name "var*k299.txt" > $RESDIR/fnames.txt"
	find $RESDIR -name "var*k299.txt" > $RESDIR/fnames.txt
	FILECNT=`find $RESDIR -name "var*k299.txt" | wc -l | cut -f1 -d' '`
	echo "Found $FILECNT files"
	
	THRESH=0
	echo "./estimateEdgeConf $RESDIR/fnames.txt $THRESH $RESDIR/edgeconf alledges"
	./estimateEdgeConf $RESDIR/fnames.txt $THRESH $RESDIR/edgeconf alledges
#done
done
#done
