DIR=~/results/networkinference/pgpm/plants/g1_g4_allsamples_initrandsearch/r4_p-5_h0.6
TOPDIR=~/results/networkinference/pgpm/plants/g1_g4_separate/
for G in g1 g2 g3 g4
do
DIR=$TOPDIR/$G/r4_p-5_h0.6
#for S in 0.2 0.4 0.6 0.8
for S in 0.1 0.2 0.3 0.4
do
mkdir $DIR/hsim${S}
echo "./reorder $DIR/coclustercnt.txt matrix  $DIR/hsim${S}/hierarchical_consensus_${S}.txt ${S} "
./reorder $DIR/coclustercnt.txt matrix  $DIR/hsim${S}/hierarchical_consensus_${S}.txt ${S} 
done
done
