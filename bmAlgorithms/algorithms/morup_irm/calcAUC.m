function [AUC,TPR,FPR]=calcAUC(West,A)
% Calculates the Area Under Curve (AUC) of the Receiver Operating
% Characteristic (ROC)
%
% Usage:
%   [AUC,TPR,FPR]=calcAUC(West,A)
% 
% Input: 
%   West:   Graph(s) of estimated probabilities of generating links in the graph for
%           entries treated as missing
%   A:      True Graph(s)
%
% Output:
%   AUC: The area under the ROC curve
%   TPR: True Positive Rate
%   FPR: False Positive Rate
%
% Written by Morten Mørup

if iscell(A)
    MAP=[];
    C=[];
    for n=1:length(A)
        [a,b,map]=find(West{n});
        MAP=[MAP; map];
        Wlogical=West{n}>0;
        [a,b,c]=find(Wlogical.*A{n}+Wlogical);
        C=[C; c-1];
    end
else
    [a,b,MAP]=find(West{1});
     Wlogical=West{1}>0;
    [a,b,C]=find(Wlogical.*A+Wlogical);
    C=C-1;
end

[val,ind]=sort(MAP(:),'ascend');
x=C(ind);
N0=sum(1-x);
N1=sum(x);
FNR=[zeros(length(x),1); 1];
TNR=[zeros(length(x),1); 1];
N_true=x(1);
N_false=1-x(1);
t=1;
for k=2:length(val)
    if val(k-1)~=val(k)
        t=t+1;
        FNR(t)=N_true/N1;
        TNR(t)=N_false/N0;
    end
    N_true=N_true+x(k);
    N_false=N_false+(1-x(k));            
end
FNR(t+1)=1;
FNR(t+2:end)=[];
TNR(t+1)=1;
TNR(t+2:end)=[];
TPR = 1-FNR;
FPR = 1-TNR;
AUC = -trapz(FPR,TPR);