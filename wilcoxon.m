
function[w] = wilcoxon(Function_name,pv,pv_name)
w=zeros(1,size(pv,1)-1);
for  i=1:size(pv,1)-1
    ansSign = signrank(pv(1,:),pv(i+1,:));     %Wilcoxon signed-rank test
    w(i)=ansSign;
end

