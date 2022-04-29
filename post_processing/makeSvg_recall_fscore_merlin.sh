

# This script will generate recall and fscore heatmaps
# usage: bash makeSvg_recall_fscore_merlin.sh
# You will  have to change file name prefixes on lines 8, 10, 24, and 26


for i in net_0.1 net_0.2 gold_standard
	do
		for j in net_0.1 net_0.2 gold_standard

			do
			if [ "$i" == "$j" ]
			then
				echo run.${i} run.${j} |awk '{printf("%s||30\t%s||30\n",$1,$2)}'
			else
				x=`../programs/validator/validate ${i}.txt ${j}.txt yes |grep -v "hit\|miss" | grep Recall | cut -f2`
				echo ${i} ${j} $x |awk '{printf("%s||30\t%s||30\t%f||%.2f\n",$1,$2,$3,$3)}'
			fi
		done
	done	| ../programs/Heatmap.awk -vC="0:(255,255,255);.50:(255,0,0)" -vL="Recall" -vD="-" > recall.svg 


for i in net_0.1 net_0.2 gold_standard
	do
			for j in net_0.1 net_0.2 gold_standard
			do
			if [ "$i" == "$j" ]
			then
				echo run.${i} run.${j} |awk '{printf("%s||30\t%s||30\n",$1,$2)}'
			else
				x=`../programs/validator/validate ${i}.txt ${j}.txt yes | grep Fscore | cut -f6`
				echo run.${i} run.${j} $x |awk '{printf("%s||30\t%s||30\t%f||%.2f\n",$1,$2,$3,$3)}'
			fi
		done
	done	| ../programs/Heatmap.awk -vC="0:(255,255,255);.50:(255,0,0)" -vL="Fscore" -vD="-" > fscore.svg 

