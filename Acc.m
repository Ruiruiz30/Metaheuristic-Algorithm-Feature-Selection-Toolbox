%_________________________________________________________________________%
%  
% Hybrid Whale Optimization Algorithm 
% with Simulated Annealing for Feature Selection 
%           By: Majdi Mafarja and Seyedali Mirjalili   
%           email: mmafarjeh@gmail.com
% 
% Main paper: M. Mafarja and S. Mirjalili                                 %
%               Hybrid Whale Optimization Algorithm                       %
%               with Simulated Annealing for Feature Selection            %
%               Neurocomputing , in press,                                %
%               DOI: https://doi.org/10.1016/j.neucom.2017.04.053         %
%                                                                         %
%  Developed in MATLAB R2014a                                             %
%                                                                         %
%  the original code of WOA is availble on                                %
%                                                                         %
%       Homepage: http://www.alimirjalili.com                             %
%                e-Mail: ali.mirjalili@gmail.com                          %
%                      
%_________________________________________________________________________%

function y=Acc(x)
%we used this function to calculate the accuracy 
global A trn vald 
 x=x>0.5;
 x=cat(2,x,zeros(size(x,1),1));
 x=logical(x);
 
if sum(x)==0
    y=[inf 
       inf];
    return;
end
mdl = fitcknn(A(trn,x),A(trn,end),'NumNeighbors',3);
c = predict(mdl,A(vald,x));
[acc] = Evaluate(A(vald,end),c);
%y=(1-SzW)*(1-cp.CorrectRate)+SzW*sum(x)/(length(x)-1);
f1=sum(x)/(length(x));
f2=1-acc(1);
y=[f1 
   f2];