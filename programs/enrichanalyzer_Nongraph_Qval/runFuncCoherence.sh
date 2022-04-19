GSIZE=2

while [ $GSIZE -le 10 ]
do
echo  "./funcCoherence ../disceq_upper_05/mbs/upper.sif ../disceq_upper_05/mbs/2n.vertex ../disceq_upper_05/mbs/slimgo/2n_slimfunc.txt ../disceq_upper_05/mbs/slimgo/func_topology_log_n10.txt 0.05 1e-6 2 $GSIZE"
echo  "./funcCoherence ../disceq_upper_05/mbs/upper.sif ../disceq_upper_05/mbs/2n.vertex ../disceq_upper_05/mbs/slimgo/2n_slimfunc.txt ../disceq_upper_05/mbs/slimgo/func_topology_log_n10.txt 0.05 1e-6 2 $GSIZE" >> run_mbs_upper_func.txt
./funcCoherence ../disceq_upper_05/mbs/upper.sif ../disceq_upper_05/mbs/2n.vertex ../disceq_upper_05/mbs/slimgo/2n_slimfunc.txt ../disceq_upper_05/mbs/slimgo/func_topology_log_n10.txt 0.05 1e-6 2 $GSIZE >> run_mbs_upper_func.txt
GSIZE=`expr $GSIZE + 1`
done
