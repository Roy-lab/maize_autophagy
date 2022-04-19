if [ $# -lt 3 ]
then
	echo "Usage: genclusterattribs dirname exprname headernames"
	exit
fi
DIR=$1
EXPR=$2
EXPRHEADERS=$3
for R in 1 2 3 4 5
do
	SUBDIR=$DIR/randinit${R}/r4_p-5_h0.6/model1/fold0
	#rm -r $SUBDIR/clusterattribs
	#mkdir -p $SUBDIR/clusterattribs
	echo "./genClusterAttrib -l $SUBDIR/goproc_chip_expr_module.txt -o $SUBDIR/clusterattribs -g $SUBDIR/goproc_details.txt -r $SUBDIR/module_reginfenr_details.txt -c $SUBDIR/module_chipenr_details.txt -d $SUBDIR/module_dnaseenr_details.txt -m $SUBDIR/modules.txt -e $EXPR -h $EXPRHEADERS"
	#./genClusterAttrib -l $SUBDIR/goproc_chip_expr_module.txt -o $SUBDIR/clusterattribs -g $SUBDIR/goproc_details.txt -r $SUBDIR/module_reginfenr_details.txt -c $SUBDIR/module_chipenr_details.txt -d $SUBDIR/module_dnaseenr_details.txt -m $SUBDIR/modules.txt -e $EXPR -h $EXPRHEADERS
	#./genClusterAttrib -l $SUBDIR/goproc_chip_expr_module.txt -o $SUBDIR/clusterattribs -g $SUBDIR/goproc_details.txt -r $SUBDIR/module_reginfenr_details.txt -c $SUBDIR/module_chipenr_details.txt -d $SUBDIR/module_koenr_details.txt -m $SUBDIR/modules.txt -e $EXPR -h $EXPRHEADERS
	#./genClusterAttrib -l $SUBDIR/goproc_chip_expr_module.txt -o $SUBDIR/clusterattribs -g $SUBDIR/goproc_details.txt -r $SUBDIR/module_reginfenr_details.txt -c $SUBDIR/module_chipenr_details.txt  -m $SUBDIR/modules.txt -e $EXPR -h $EXPRHEADERS
done
