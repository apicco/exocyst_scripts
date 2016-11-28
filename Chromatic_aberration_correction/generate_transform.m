!grep -v [a-z] RFP_0 | sed 's/,//g' > tmpW1
!grep -v [a-z] GFP_0 | sed 's/,//g' > tmpW2
rawW1=load('tmpW1'); 
rawW2=load('tmpW2'); 
!rm tmpW[12]

cutoff=2
W1=[];
W2=[];
d=[];
N=12
for i = 1:length(rawW1(:,1))
	dtmp=sqrt((rawW1(i,1)-rawW2(:,1)).^2+(rawW1(i,2)-rawW2(:,2)).^2);
	if (length(find(dtmp==min(dtmp)))==1 && min(dtmp)<cutoff)
		W1=[W1; rawW1(i,:)];
		W2=[W2; rawW2(find(dtmp==min(dtmp)),:)];
		d=[d;min(dtmp)];
	end
end

save('W1_ref.data','-ascii','W1')
save('W2_ref.data','-ascii','W2')

warp=cp2tform(W2(:,1:2),W1(:,1:2),'lwm',N);
save('transform','warp')
