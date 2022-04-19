if [ $# -lt 2 ]
then
	echo "./batchAssess INDIR OUTDIR"
	exit
fi


for CONS in 0.0 0.2 0.4 0.6 0.8
do
	for MOECNT in 10 15 20 25
	do
		INDIR=$1/cons$CONS/moe/k$MOECNT
		echo "find $INDIR -name \"clusterassign*\" > scratch/filelist_cons${CONS}_k${MOECNT}.txt"
		find $INDIR -name "clusterassign*" > scratch/filelist_cons${CONS}_k${MOECNT}.txt
		ADJFNAME=$2/adj_cons${CONS}_k${MOECNT}.txt
		echo "bsub -o scratch/bsub$RUN.txt -q long ./assessClusterStab scratch/filelist_cons${CONS}_k${MOECNT}.txt $ADJFNAME"
		#bsub -o scratch/bsub$RUN.txt -q long ./assessClusterStab scratch/filelist_cons${CONS}_k${MOECNT}.txt $ADJFNAME
	done
done
