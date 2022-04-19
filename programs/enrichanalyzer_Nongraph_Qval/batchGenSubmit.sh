if [ $# != 1 ]
then
        echo "Incorrect number of arguments $#"
        echo "Usage: genSubmitFile submitfilename "
        exit
fi
SUBMITFNAME=$1
echo "universe = vanilla" >  $SUBMITFNAME
echo "executable = wrapper_enrichanalyzer.sh" >> $SUBMITFNAME

#for K in 7 9 11 13 15 17 19 21 23 25
for K in 5 7 9 11 13 15 17 19 21 23 25
do
for RANDINIT in 1 2 3 4 5 6 7 8 9 10
do
OUTDIR=scratch_cichlids/k${K}/run${RANDINIT}/randinit${RANDINIT}
RUNDIR=scratch_cichlids/k${K}/run${RANDINIT}
ACTUALRESULTDIR=../../../results/compfuncgen/cross_speciescluster/cichlids/k${K}/randinit${RANDINIT}
RESULTDIR=../../../../../../results/compfuncgen/cross_speciescluster/cichlids/k${K}/randinit${RANDINIT}
if [  -f $ACTUALRESULTDIR/Ga_clusterassign.txt  ]
then
echo "mkdir -p $OUTDIR"
mkdir -p $OUTDIR
for SPECIES in  Ga Mz Nb On Pn Ab Anc1 Anc2 Anc3 Anc4 Anc5
do
echo "error = jobs_k${K}_r${RANDINIT}.err" >> $SUBMITFNAME
echo "log = jobs_k${K}_r${RANDINIT}.log" >> $SUBMITFNAME
echo "output = jobs_k${K}_r${RANDINIT}.out" >> $SUBMITFNAME
echo "initialdir = $RUNDIR" >> $SUBMITFNAME
echo "match_list_length = 5" >> $SUBMITFNAME
echo "requirements = (TARGET.Name =!= LastMatchName1) && ( OpSysAndVer=="SL6" || OpSysAndVer=="RedHat6" || IsRHEL6=="True")" >> $SUBMITFNAME
echo "Rank=kflops" >> $SUBMITFNAME
echo "should_transfer_files = YES" >> $SUBMITFNAME
echo "+WantRHEL6Job = true" >> $SUBMITFNAME 
echo "+WantFlocking = true" >> $SUBMITFNAME 
echo "+group=\"WID\"" >> $SUBMITFNAME
echo "when_to_transfer_output = ON_EXIT" >> $SUBMITFNAME
echo "Notification = never" >> $SUBMITFNAME
echo "arguments = \"$SPECIES \"" >> $SUBMITFNAME  
echo "transfer_input_files = ../../../enrichAnalyzer,../../../../../../data/cichilddata/geneontology/species_expanded/on_go_regnet.txt,$RESULTDIR/${SPECIES}_clusterassign.txt,$RESULTDIR/${SPECIES}_genesets.txt,../../../sharedlib.tgz" >> $SUBMITFNAME
echo "transfer_output_files = ${SPECIES}_goenr_details.txt" >> $SUBMITFNAME
echo "transfer_output_remaps = \"${SPECIES}_goenr_details.txt= $RESULTDIR/${SPECIES}_goenr_details.txt\"" >> $SUBMITFNAME
echo "queue" >> $SUBMITFNAME
echo"" >> $SUBMITFNAME
done
fi
done
done
