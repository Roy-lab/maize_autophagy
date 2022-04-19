
if [ $# -ne 1 ] 
then
	echo "Usage batchEstimateConf outputdir"
	#exit
fi
DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/stressdata/ilan_carbonstress_pp/
DIR=/seq/compbio/sroy/networklearning/softclustering/clustering/ilandata_run/carbondata/
K=5
while [ $K -le 5 ]
do	
	for LEAF in  0.8
	do
	for RITER in 1 2 3 4 5 6 7 8 9 10
	do
	#RESDIR=$DIR/crossspecies_k${K}/randinit${RITER}/lf${LEAF}_nlf${LEAF}
	RESDIR=$DIR/clusteroutput_k${K}/randinit${RITER}/iter2
	#RESDIR=$DIR/crossspecies_k${K}/randinit${RITER}/
	#for SPECIES in `cat $SPECLIST`
	for SPECIES in Scer Cgla Scas Klac Calb
	do
		SUBSTR=Anc
		LEN=`expr match $SPECIES $SUBSTR`
		if [ $LEN -gt 0 ]
		then
			echo "Continuing at $SPECIES $LEN"
			continue
		fi
		OUTSUFF=$RESDIR/${SPECIES}_goprocenr.txt
	done
done
