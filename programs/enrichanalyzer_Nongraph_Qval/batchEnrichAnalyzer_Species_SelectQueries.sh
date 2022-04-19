TOPDIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/carbonclusters_dup2c/allclade4/crossspecies_k5/randinit4/lf0.8_nlf0.8/
#TOPDIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/carbon_mergedassign_wkp3/allclade4/crossspecies_k5/randinit1/lf0.8_nlf0.8/
#DIR=/ahg/regev/users/sroy/compfuncgen/fordawn/selectedqueries/formatedlists_new/speciesspecific_lists
#TOPDIR=/ahg/regev/users/sroy/compfuncgen/fordawn/sfp1/
MOTIFDIR=/seq/compbio/sroy/networklearning/data/motifs_t0.6
MOTIFDIR=/ahg/regev/users/sroy/compfuncgen/data/Motifs_MacIsaac06/
GODIR=/ahg/regev/users/sroy/compfuncgen/data/speciespec_go/
SPECLIST=../scripts/speclist4_wkp.txt
SPECLIST=../scripts/speclist4.txt
GOFILE=/seq/compbio/sroy/networklearning/data/go/goproc
#SPECLIST=../scripts/speclist_heat_8spec.txt
#LISTDIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/carbonclusters_dup2c/allclade4/crossspecies_k5/randinit4/lf0.8_nlf0.8
#TOPDIR=/seq/compbio/sroy/networklearning/fgraph/results/gmm_crosspecies/stressdata/ilan_heatstress_8spec_calb42_cgla37_pp_postem_2/crossspecies_k5/randinit2/lf0.8_nlf0.8/
LISTDIR=/ahg/regev/users/sroy/compfuncgen/fordawn/
LISTDIR=/ahg/regev/users/sroy/compfuncgen/data/sfp1/
#EXPRDIR=/ahg/regev/users/sroy/compfuncgen/data/carbon/species_specific
#for SUBSTR in lca_prepost_pattens/lca_diverged lca_prepost_patterns/diverged lca_prepost_patterns/left_conserved lca_prepost_patterns/right_conserved
#for SUBSTR in lca_schizo_nonschizo/left_conserved lca_schizo_nonschizo/right_conserved
#for SUBSTR in gene_sets
#for SUBSTR in cladespecific_genesets avivselect_genesets
#for SUBSTR in querylist_manual_novel3
for SUBSTR in genesets
#for SUBSTR in diff_exp_genesets
do

DIR=$LISTDIR/$SUBSTR/
#DIR=$TOPDIR/speciesspecific_lists
ITER=1
#for SPECIES in `cat $SPECLIST`
for SPECIES in   Cgla Scas  Klac Spom
#for SPECIES in  Scer Cgla1 Cgla2 Scas1 Scas2  Klac Spom
do
	SUBSTR=Anc
	LEN=`expr match $SPECIES $SUBSTR`
	if [ $LEN -gt 0 ]
	then
		echo "Continuing at $SPECIES $LEN"
		continue
	fi
	MOTIFFILE=$MOTIFDIR/${SPECIES}_regnet.txt
	#GOFILE=$GODIR/${SPECIES}_goproc
	#OUTSUFF=$DIR/${SPECIES}_macissac_motifenr
	OUTSUFF=$DIR/${SPECIES}_goenr_Plus2_withhits
	#NWFILE=$DIR/${SPECIES}_genesets.txt
	NWFILE=$DIR/${SPECIES}_genesets_Plus2.txt
	#cat `find  $DIR/$SPECIES/*.txt` | sort -u > $DIR/$SPECIES/bg.list
	#GENELIST=$TOPDIR/${SPECIES}_clusterassign.txt
	#GENELIST=$TOPDIR/${SPECIES}_speciesspecnames_clusterassign.txt
	#GENELIST=$DIR/${SPECIES}.txt
	GENELIST=$DIR/${SPECIES}_clusterassign_Plus2.txt
	SPLOWER=`echo $SPECIES | awk '{printf("%s\n",tolower($1))}'`
	#GENELIST=$EXPRDIR/$SPLOWER/nointerpolate_0.50miss.geneexp
	#GENELIST=$DIR/${SPECIES}/bg.list
	#echo "bsub -o scratch/selectqueries_$ITER.txt ./enrichAnalyzer $NWFILE $GENELIST $MOTIFFILE 0.05 $OUTSUFF persg"
	#bsub -q hour -o scratch/selectqueries_$ITER.txt ./enrichAnalyzer $NWFILE $GENELIST $MOTIFFILE 0.05 $OUTSUFF persg
	echo "bsub -o scratch/selectqueries_$ITER.txt ./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $OUTSUFF persg"
	bsub -q hour -o scratch/selectqueries_$ITER.txt ./enrichAnalyzer $NWFILE $GENELIST $GOFILE 0.05 $OUTSUFF persg
	ITER=`expr $ITER + 1`
done
done
