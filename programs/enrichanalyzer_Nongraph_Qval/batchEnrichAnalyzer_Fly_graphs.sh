IMAGOFILE=/seq/compbio/sroy/modencode/fly/flynet/data/imago/fbgn_imago_slimterm.txt
#RANDFILE=/seq/compbio/sroy/networklearning/data/randgo/randgoproc_s3-44_n1000_k20.txt
GENELIST=/seq/compbio/sroy/modencode/fly/flynet/overlaps_v3/nodes/Cr_Ch_Mo_Im_Rseq_tftgt.genelist
#GENELIST=/seq/compbio/sroy/modencode/fly/exprdata/suedata/Cr_Ch_Mo_Im_SueTC.genelist

OUTPUTDIR=/seq/compbio/sroy/modencode/fly/analysis/priorbased_Cr_Ch_Mo_Im_Rseq
RUN=1
#for DIRSUFF in normaffinity/cons0.0 normaffinity/cons0.2 normaffinity/cons0.4 normaffinity/cons0.6 normaffinity/cons0.8 
#for DATA in motif
for DATA in motif chip chromatin motif_chip motif_chip_chromatin
#for DATA in motif_symm chip_symm chromatin_symm motif_chip_symm motif_chip_chromatin_symm
do
for BETA1 in 3 4
do
	for BETA2 in 3 4
	do
	RESDIR=$OUTPUTDIR/$DATA/beta${BETA1}_${BETA2}/model1/fold0
	for K in 8 10 12 14 16 28
	do
	NWFILE=$RESDIR/genesets${K}.txt
	OUTSUFF=$RESDIR/imago_unsign_k${K}
	SCRATCHFILE=scratch/fly_${DATA}_symm_imago_${RUN}.txt
	if [ -f $SCRATCHFILE ]
	then
		echo "$SCRATCHFILE"
		rm $SCRATCHFILE
	fi
	echo "bsub -q short -o $SCRATCHFILE ./enrichAnalyzer $NWFILE $GENELIST $IMAGOFILE 0.05 $OUTSUFF persg"
	bsub -q short -o  $SCRATCHFILE ./enrichAnalyzer $NWFILE $GENELIST $IMAGOFILE 0.05 $OUTSUFF persg
	RUN=`expr $RUN + 1`
	done
	done
	done
done
