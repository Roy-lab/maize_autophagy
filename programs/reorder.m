function reorder(input_mat)

	%a=importdata('/mnt/dv/wid/projects5/Roy-singlecell/sgm_work/otegui/merlin54/output/processed/merged3mats_hvg2folds.txt')
	a=importdata(imput_mat);
	h1=a.textdata(1,1);
	header=a.textdata(1,2:end);
	genenames=a.textdata(2:end,1);
	d=a.data;
	fprintf('dist\n');
	dist=pdist(d');
	fprintf('tree\n');
	tree=linkage(dist,'average');
	fprintf('reorder\n');
	ord=optimalleaforder(tree,dist);
	fprintf('reorder\n');
	dd=d(:,ord);
	fprintf('write\n');
	hh=header(ord);
	newheader=[h1,hh];
	saveTab('reordered_exp.txt',newheader,genenames,dd)
