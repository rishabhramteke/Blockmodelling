function [distance] = varOfInfoUnlabel(mMatchingWeights, vPosWeight1, vPosWeight2)
%
% Computes the variation of information for external cluster comparison.
%
% vMembership1 - vector of membership (each element is 1..K, where K is the
% number of clusters in clustering set 1)
% vMembership2 - vector of membership (each element is 1..K', where K' is the
% number of clusters in clustering set 2)
%


vDistrib1 = vPosWeight1;
vDistrib2 = vPosWeight2;
mInterDistrib = zeros(size(vPosWeight1,1), size(vPosWeight2,1));
for r = 1 : size(vPosWeight1,1)
    for c = 1 : size(vPosWeight2,1)
        %mInterDistrib(r,c) = vPosWeight1(r) * mMatchingWeights((r-1)*size(vPosWeight2,1) + c) * vPosWeight2(c);
        mInterDistrib(r,c) = mMatchingWeights((r-1)*size(vPosWeight2,1) + c);
    end
end


entropyVal1 = myEntropy(vDistrib1);
entropyVal2 = myEntropy(vDistrib2);

mutualVal = 0;
for i = 1 : size(vPosWeight1,1)
    for j = 1 : size(vPosWeight2,1)
        mutualVal = mutualVal + mInterDistrib(i,j) * myLog(mInterDistrib(i,j) / (vDistrib1(i) * vDistrib2(j)));
    end
end

distance = entropyVal1 + entropyVal2 - 2 * mutualVal;

end



function [logVal] = myLog(val)

if val <= 0
    logVal = 0;
else
    logVal = log2(val);
end

end



function [value] = myEntropy(vDistribution)
%
% compute entropy
%


value = 0;
for p = 1 : size(vDistribution,1)
    if vDistribution(p) > 0
        value = value + vDistribution(p) * log2(vDistribution(p));
    end
end
value = -value;

end