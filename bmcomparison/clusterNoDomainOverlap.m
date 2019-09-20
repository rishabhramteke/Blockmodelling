function clusterNoDomainOverlap(cSnapshots, cMemberships, bPosDef, sMatchType, sumWeight)
%
% Cluster the input blockmodels
%

mPosDis = zeros(size(cSnapshots,1), size(cSnapshots,1));
mDenDis = zeros(size(cSnapshots,1), size(cSnapshots,1));

% compute the pairwise distance measures
for i = 1 : size(cSnapshots,1)
    for j = i+1 : size(cSnapshots,1)
        [posDis, densDis, mW] = compareNoDomainOverlap(cMemberships{i}, cSnapshots{i}, cMemberships{j}, cSnapshots{j}, bPosDef, sMatchType);
        mPosDis(i,j) = posDis;
        mDenDis(i,j) = densDis;
    end
end

% use weighted sum to combine the two distances
mCombDis = sumWeight * mPosDis + (1-sumWeight) * mDenDis;

% perform the clustering



end