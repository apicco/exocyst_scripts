
load('transform.mat')
W1files=dir('RFP_0');
W2files=dir('GFP_0');
for i=1:length(W1files)
	rawW1=load(W1files(i).name); 
	rawW2=load(W2files(i).name); 
	
	newW1=tforminv(warp,rawW1(:,1:2));
	outputW1=[newW1,rawW1(:,3:length(rawW1(1,:)))];
	save(strcat(W1files(i).name(1:(length(W1files(i).name))),'_warped'),'-ascii','outputW1')
end

exit
