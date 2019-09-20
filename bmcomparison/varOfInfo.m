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
% map each symbol to a continuous range

[r1,c1] = size(vSymbols1);
numRVValues1 = r1;
if (c1 > r1)
    numRVValues1 = c1;
end
[r2,c2] = size(vSymbols2);
numRVValues2 = r2;
if (c2 > r2)
    numRVValues2 = c2;
end

vSymbolMap1 = zeros(1, max(vSymbols1));
for i = 1 : numRVValues1
    vSymbolMap1(vSymbols1(i)) = i;
end

vSymbolMap2 = zeros(1, max(vSymbols2));
for i = 1 : numRVValues2
    vSymbolMap2(vSymbols2(i)) = i;
end

% compute the distributions
vDistrib1 = distribution(vMembership1, vSymbols1, vSymbolMap1);
vDistrib2 = distribution(vMembership2, vSymbols2, vSymbolMap2);
mInterDistrib = intersectDistrib(vMembership1, vMembership2, vSymbols1, vSymbols2, vSymbolMap1, vSymbolMap2);

entropyVal1 = myEntropy(vDistrib1);
entropyVal2 = myEntropy(vDistrib2);

mutualVal = 0;
dim = 1;
if (size(vMembership1,2) > size(vMembership1,1))
    dim = 2;
end
for i = 1 : size(vSymbols1,dim)
    for j = 1 : size(vSymbols2,dim) 
        % mylog returns negative values, so we adding negative values, so
        % we need to negate it first
        mutualVal = mutualVal + mInterDistrib(i,j) * myLog(mInterDistrib(i,j) / (vDistrib1(i) * vDistrib2(j)));
    end
end

distance = entropyVal1 + entropyVal2 - 2 * mutualVal;

end


function [mInterDistrib] = intersectDistrib(vMembership1, vMembership2, vSymbols1, vSymbols2,...
    vSymbolMap1, vSymbolMap2)
%
% Compute the intersecting distribution for two cluster membership vectors.
%

dim = 1;
if (size(vMembership1,2) > size(vMembership1,1))
    dim = 2;
end

assert(size(vMembership1,dim) == size(vMembership2,dim));

% use the bigger dimension
[r1,c1] = size(vSymbols1);
numRVValues1 = r1;
if (c1 > r1)
    numRVValues1 = c1;
end
[r2,c2] = size(vSymbols2);
numRVValues2 = r2;
if (c2 > r2)
    numRVValues2 = c2;
end


mInterDistrib = zeros(numRVValues1, numRVValues2);

cClusIndices1 = cell(1, numRVValues1);
cClusIndices2 = cell(1, numRVValues2);

for i = 1 : numRVValues1
    cClusIndices1{i} = find(vMembership1 == vSymbols1(i));
end

for j = 1 : numRVValues2
    cClusIndices2{j} = find(vMembership2 == vSymbols2(j));
end


for i = 1 : numRVValues1
    for j = 1 : numRVValues2
        vOverlap = intersect(cClusIndices1{i}, cClusIndices2{j});
        overlapSize = size(intersect(cClusIndices1{i}, cClusIndices2{j}), dim);
        if isempty(vOverlap)
            overlapSize = 0;
        end
        mInterDistrib(i,j) = overlapSize / size(vMembership1,dim);
    end
end
 

end



function [vDistribution] = distribution(vElems, vSymbols, vSymbolMap)
%
% vDistribution - row vector (1 x n) 
%


% figure out with orientation is vElems and vSymbols
[r,c] = size(vSymbols);
numRVValues = r;
if (c > r)
    numRVValues = c;
end
vDistribution = zeros(1,numRVValues);





% check which dimension we are iterating over
% it should be the larger dimension
dim = 1;
if size(vElems,2) > size(vElems,1)
    dim = 2;
end

% compute distribution
for i = 1 : size(vElems,dim)
   vDistribution(vSymbolMap(vElems(i))) = vDistribution(vSymbolMap(vElems(i))) + 1;
end

% normalise to frequency
for s = 1 : size(vDistribution,2)
   vDistribution(s) = vDistribution(s) / size(vElems,dim); 
end

end


% function [logVal] = myLog(val)
% 
% if val == 0
%     logVal = 0;
% else
%     logVal = log2(val);
% end
% 
% end


% 
% function [value] = myEntropy(vDistribution)
% %
% % compute entropy
% %
% 
% 
% value = 0;
% for p = 1 : size(vDistribution,2)
%     if vDistribution(p) > 0
%         value = value + vDistribution(p) * log2(vDistribution(p));
%     end
% end
% 
% value = -value;
% 
% end