

for i in $1; 
do 
	mkdir -p heatmaps_out/$i
	for j in `ls heatmaps_in/$i |grep _attrib.txt |sed 's/_attrib\.txt//g'`; 
	do 
		echo $i.$j; 
		cat heatmaps_in/$i/${j}_attrib.txt heatmaps_in/$i/${j}_regulators.txt | ./Heatmap.awk -vC="0:(255,255,255);1:(0,0,255);2:(0,255,255);3:(0,100,100);4:(100,0,100);5:(80,0,255) -2:(0,0,255);0:(255,255,255);2:(255,0,0)" -vD=" " -vFontSize=7 -vStrokeC="-" -vStokeW=0.5 -vStrokeSC="black"  -vL="regulators expression"  > heatmaps_out/$i/$j.svg; 
		convert heatmaps_out/$i/$j.svg heatmaps_out/$i/$j.png
	done; 
	#cat heatmaps_in/$i/ModuleAvg.txt | ./Heatmap.awk -vC="-2:(0,0,255);0:(255,255,255);2:(255,0,0)" -vD=" " -vFontSize=7 -vStrokeC="-" -vStokeW=0.5 -vStrokeSC="black"  -vL="avgExpression"  > heatmaps_out/$i/ModuleAvg.svg; 
	cat heatmaps_in/$i/ModuleAvg.txt | ./Heatmap.awk -vC="-1:(0,0,255);0:(255,255,255);1:(255,0,0)" -vD=" " -vFontSize=7 -vStrokeC="-" -vStokeW=0.5 -vStrokeSC="black"  -vL="avgExpression"  > heatmaps_out/$i/ModuleAvg.svg; 
	convert heatmaps_out/$i/ModuleAvg.svg heatmaps_out/$i/ModuleAvg.png
done
