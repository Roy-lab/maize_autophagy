if [ $# -lt 1 ]
then
	echo "./batchOptimal OUTDIR"
	exit
fi

OPTIMAL=/ahg/regev/users/sroy/thirdparty/Optimal/arrange

for CONS in 0.0 0.2 0.4 0.6 0.8
do
	for MOECNT in 10 15 20 25
	do
		ADJFNAME=$1/adj_cons${CONS}_k${MOECNT}
		echo "bsub -o scratch/bsub$RUN.txt -q long $OPTIMAL $ADJFNAME"
		bsub -o scratch/bsub$RUN.txt -q long $OPTIMAL $ADJFNAME
		RUN=`expr $RUN  + 1`
	done
done
