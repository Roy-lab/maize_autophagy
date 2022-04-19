if [ $# -ne 1 ] 
then
	echo "Usage batchEstimateConf outputdir"
	exit
fi
MOTIFDIR=~/data/reprogramming/rupadata/binding
K=15
DIR=$1
while [ $K -le 15 ]
do	
	for RITER in 1 2 3 4 5 
	do
	RESDIR=$DIR/k${K}/randinit${RITER}/
	OUTDIR=$RESDIR/
	for SPECIES in ips pips mef
	do
		SUBSTR=Anc
		LEN=`expr match $SPECIES $SUBSTR`
		if [ $LEN -gt 0 ]
		then
			echo "Continuing at $SPECIES $LEN"
			continue
		fi
		NWFILE=$RESDIR/${SPECIES}_genesets.txt
		MOTIFFILE=$MOTIFDIR/${SPECIES}_regnet.txt
		GENELIST=$RESDIR/${SPECIES}_clusterassign.txt
		OUTSUFF=$RESDIR/${SPECIES}_motifenr
		echo "./enrichAnalyzer $NWFILE $GENELIST $MOTIFFILE 0.05 $OUTSUFF persg"
		./enrichAnalyzer $NWFILE $GENELIST $MOTIFFILE 0.05 $OUTSUFF persg
	#done
	done
	done
	K=`expr $K + 2`
done

