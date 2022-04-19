

DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/carbon_mergedassign_wkp3/allclade4/crossspecies_k5/randinit1/lf0.8_nlf0.8
SPECLIST=../scripts/speclist4_wkp.txt
SUBDIR=$DIR/species_specific_termlist
mkdir -p $SUBDIR

for SPECIES in `cat $SPECLIST`
do
	SUBSTR=Anc
	LEN=`expr match $SPECIES $SUBSTR`
	if [ $LEN -gt 0 ]
	then
		echo "Continuing at $SPECIES $LEN"
		continue
	fi
	FNAME=$DIR/${SPECIES}_subgenesets_enrgoproc_${SPECIES}.txt
	echo "cut -f1 $FNAME | cut -f2 -d':' | sort -u > $SUBDIR/${SPECIES}_termlist.txt"
	cut -f1 $FNAME | cut -f2 -d':' | sort -u > $SUBDIR/${SPECIES}_termlist.txt
done
