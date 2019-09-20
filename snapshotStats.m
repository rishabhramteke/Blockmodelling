function [vEdgeOverLap, vEdgePercPred,vEdgePercSucc, vEdgePercUnion, vNonEdgePercUnion ] = snapshotStats(cmSnapshots, cPartitions, aPartInfo)
%
% Computes various statistics to do with the snapshots and discovered
% clusters.
%
% cmSnapshots - a cell vector of 2d matrices of the snapshots.
% cPartitions - a cell of vectors of the partitions (indices into the
% snapshots).  Each row is the partitions of one subsequence
% aPartInfo - an array of partition info, each row is the information of
% one subsequence
%

% make sure the sizes of the snapshots are the same.
assert(size(aPartInfo,1) == size(cPartitions,1));
% make sure the length of the subsequences is the same as the number of
% snapshots
assert(size(cmSnapshots,2) == sum(aPartInfo(end,2)));

vEdgeOverLap = zeros(1, size(cmSnapshots,2)-1);
vEdgePercPred = zeros(1, size(cmSnapshots,2)-1);
vEdgePercSucc = zeros(1, size(cmSnapshots,2)-1);
vEdgePercUnion = zeros(1, size(cmSnapshots,2)-1);
vNonEdgePercUnion = zeros(1, size(cmSnapshots,2)-1);

for i = 1 : size(cmSnapshots,2) - 1
    vEdgeOverLap(i) = edgeOverlap(cmSnapshots{i}, cmSnapshots{i+1});
    % number of edges in i
    nzEdges1 = sum(nonzeros(cmSnapshots{i}));
    vEdgePercPred(i) = vEdgeOverLap(i) / nzEdges1;
    % number of egges in i+1
    nzEdges2 = sum(nonzeros(cmSnapshots{i+1}));
    vEdgePercSucc(i) = vEdgeOverLap(i) / nzEdges2;
    % number of edges in both
    nzEdgesUnion = edgeUnionNum(cmSnapshots{i}, cmSnapshots{i+1});
    vEdgePercUnion(i) = vEdgeOverLap(i) / nzEdgesUnion;
    % non edges
    zEdgesUnion = size(union(find(cmSnapshots{i} == 0), find(cmSnapshots{i+1} == 0)));
    vNonEdgePercUnionA
    zEdgesUnion = size(union(find(cmSnapshots{i} == 0), find(cmSnapshots{i+1} == 0)));
    vNonEdgePercUnionA(i) = nonEdgeOverlap(cmSnapshots{i}, cmSnapshots{i+1}) / zEdgesUnion;
end


end % end of function


function [overlapVal] = edgeOverlap(mAdj1, mAdj2)
%
% Compute the overlap between mAdj1 and mAdj2.
%

overlapVal = sum(nonzeros(mAdj1 .* mAdj2));

end % end of function


function [overlapVal] = nonEdgeOverlap(mAdj1, mAdj2)
%
% Compute the overlap between mAdj1 and mAdj2.
%

%overlapVal = sum(nonzeros(1 - xor(mAdj1,mAdj2))) - edgeOverlap(mAdj1, mAdj2)
overlapVal = size(intersect(find(mAdj1 == 0), find(mAdj2 == 0)));

end % end of function


function [unionVal] = edgeUnionNum(mAdj1, mAdj2)
%
% Compute the overlap between mAdj1 and mAdj2.
%

unionVal = sum(nonzeros(mAdj1 + mAdj2));

end % end of function