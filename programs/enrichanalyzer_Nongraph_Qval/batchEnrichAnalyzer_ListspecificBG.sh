
DIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/carbon_mergedassign_wkp3/allclade4/crossspecies_k5/randinit1/lf0.8_nlf0.8
MOTIFDIR=/seq/compbio/sroy/networklearning/data/motifs_t0.6
BGLISTDIR=/ahg/regev/users/sroy/compfuncgen/data/speciespec_go/
SPECLIST=../scripts/speclist4_wkp.txt

f#or SPECIES in `cat $SPECLIST`
for SPECIES in Scer
do
FNAME=$DIR/${SPECIES}_subgenesets_enrgoproc_${SPECIES}.txt
TERMLIST=$DIR/species_specific_termlist/${SPECIES}_termlist.txt
mkdir $SUBDIR
for TERM in `cat ${TERMLIST}`
do
	echo $TERM
	export TERM
	echo "awk '{m=match($1,ENVIRON["TERM"]); if(m>0) printf("%s\n",$0)}' $FNAME > scratch_listspecific_bg/${SPECIES}_${TERM}_query.txt"
	awk '{m=match($1,ENVIRON["TERM"]); if(m>0) printf("%s\n",$0)}' $FNAME > scratch_listspecific_bg/${SPECIES}_${TERM}_query.txt
	NWFILE=scratch_listspecific_bg/${SPECIES}_${TERM}_query.txt
	BGLIST=$BGLISTDIR/${SPECIES}_lists/${TERM}.txt
	OUTPUTSUFF=listspecific_bg_enrichment/${SPECIES}_${TERM}_motifenr
	echo "bsub -q hour -o scratch_listspecific_bg/run_${SPECIES}_${TERM}.txt ./enrichAnalyzer $NWFILE $BGLIST  $MOTIFFILE 0.05 $OUTSUFF persg"
	bsub -q hour -o scratch_listspecific_bg/run_${SPECIES}_${TERM}.txt ./enrichAnalyzer $NWFILE $BGLIST  $MOTIFFILE 0.05 $OUTSUFF persg
done
done
