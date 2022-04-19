
function zeromean(input_mat)

	a=importdata(imput_mat);
	h1=a.textdata(1,1);
	header=a.textdata(1,2:end);
	genenames=a.textdata(2:end,1);
	d=a.data;
	m=mean(d,2);
	m=repmat(m,1,size(d,2));
	z=d-m;
	newheader=[h1,header];
	saveTab('exp_reordered_zeromean.txt',newheader,genenames,z)


quit

