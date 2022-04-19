IMAGOFILE=/seq/compbio/sroy/modencode/fly/flynet/data/imago/fbgn_imago_slimterm.txt
IMAGOFILE=/seq/compbio/sroy/modencode/fly/flynet/data/go/flygoproc
MOTIFILE=/seq/compbio/sroy/modencode/fly/motifdata/fly_7mer_0.6_regnet.txt
#RANDFILE=/seq/compbio/sroy/networklearning/data/randgo/randgoproc_s3-44_n1000_k20.txt
#GENELIST=/seq/compbio/sroy/modencode/fly/exprdata/kwrpkm_0.80miss_high1_std1.genelist
#GENELIST=/seq/compbio/sroy/modencode/fly/exprdata/v2/imago_female.genelist

OUTPUTDIR=/seq/compbio/sroy/modencode/fly/analysis/gmmclusters/
RUN=1
#for DIRSUFF in normaffinity/cons0.0 normaffinity/cons0.2 normaffinity/cons0.4 normaffinity/cons0.6 normaffinity/cons0.8 
#for DIRSUFF in normaffinity/cons0.0
#for DIRSUFF in  KW/k34 SC_withgamma/k34
for DIRSUFF in  SC_withgamma/k34
do
ITER=1
while [ $ITER -le 10 ]
do
	RESDIR=$OUTPUTDIR/$DIRSUFF/randinit${ITER}
	NWFILE=$RESDIR/genesets.txt
	GENELIST=$RESDIR/clusterassign.txt
	OUTSUFF=$RESDIR/motifenr
	echo "bsub -o scratch/fly_cons_go_${ITER}_${RUN}.txt ./enrichAnalyzer $NWFILE $GENELIST $IMAGOFILE 0.05 $OUTSUFF persg"
	SCRATCHFILE=scratch/fly_cons_motif_${ITER}_${RUN}.txt
	if [ -f $SCRATCHFILE ]
	then
		rm $SCRATCHFILE
	fi
	#bsub -o $SCRATCHFILE ./enrichAnalyzer $NWFILE $GENELIST $IMAGOFILE 0.05 $OUTSUFF persg
	bsub -o $SCRATCHFILE ./enrichAnalyzer $NWFILE $GENELIST $MOTIFILE 0.05 $OUTSUFF persg
	ITER=`expr $ITER + 1`
	RUN=`expr $RUN + 1`
done
done
