function saveTab(outname,header,names,q)

f=fopen(outname,'W');
fprintf(f,'%s',header{1});
for i=2:length(header)
	fprintf(f,'\t%s',header{i});
end
fprintf(f,'\n');
for i=1:size(q,1)
	fprintf(f,'%s',names{i});
	fprintf(f,'\t%f',q(i,:));
	%for j=1:size(q,2)
	%	fprintf(f,'\t%f',q(i,j));
	%end
	fprintf(f,'\n');
end
fclose(f);

