function [distance] = varOfInfo(vMembership1, vMembership2)
%
% Computes the variation of information for external cluster comparison.
%
% vMembership1 - vector of membership (each element is 1..K, where K is the
% number of clusters in clustering set 1)
% vMembership2 - vector of membership (each element is 1..K', where K' is the
% number of clusters in clustering set 2)
%

vSymbols1 = unique(vMembership1);
vSymbols2 = unique(vMembership2);

vDistrib1 = distribution(vMembership1)
vDistrib2 = distribution(vMembership2)
mInterDistrib = intersectDistrib(vMembership1, vMembership2)

entropyVal1 = myEntropy(vDistrib1)
entropyVal2 = myEntropy(vDistrib2)

mutualVal = 0;
for i = 1 : size(vSymbols1,2)
    for j = 1 : size(vSymbols2,2) 
        mutualVal = mutualVal + mInterDistrib(i,j) * myLog(mInterDistrib(i,j) / (vDistrib1(i) * vDistrib2(j)))
    end
end

distance = entropyVal1 + entropyVal2 - 2 * mutualVal;

end


function [mInterDistrib] = intersectDistrib(vMembership1, vMembership2)
%
% Compute the intersecting distribution for two cluster membership vectors.
%

assert(size(vMembership1,2) == size(vMembership2,2));


vSymbols1 = unique(vMembership1);
vSymbols2 = unique(vMembership2);

mInterDistrib = zeros(size(vSymbols1,2), size(vSymbols2,2));

cClusIndices1 = cell(1, size(vSymbols1,2));
cClusIndices2 = cell(1, size(vSymbols2,2));

for i = 1 : size(vSymbols1,2)
    cClusIndices1{i} = find(vMembership1 == i);
end

for i = 1 : size(vSymbols2,2)
    cClusIndices2{i} = find(vMembership2 == i);
end


for i = 1 : size(vSymbols1,2)
    for j = 1 : size(vSymbols2,2) 
        mInterDistrib(i,j) = size(intersect(cClusIndices1{i}, cClusIndices2{j}), 2) / size(vMembership1,2);
    end
end


end



function [vDistribution] = distribution(vElems)

vSymbols = unique(vElems);
vDistribution = zeros(1,size(vSymbols,2));

% compute distribution
for i = 1 : size(vElems,2)
   vDistribution(vElems(i)) = vDistribution(vElems(i)) + 1;
end

% normalise to frequency
for s = 1 : size(vDistribution,2)
   vDistribution(s) = vDistribution(s) / size(vElems,2); 
end

end


function [logVal] = myLog(val)

if val == 0
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
for p = 1 : size(vDistribution,2)
    if vDistribution(p) > 0
        value = value + vDistribution(p) * log2(vDistribution(p));
    end
end

value = -value;

end