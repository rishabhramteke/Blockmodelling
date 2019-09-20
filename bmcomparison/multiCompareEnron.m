function [mDistances] = multiCompareEnron(cSnapshots, cPartitions, aPartInfo)
%
% Compares the conseqative blockmodels in the sequence of cPartitions.
%

assert(size(cSnapshots,2) == size(cPartitions,1));
assert(size(aPartInfo, 1) == size(cPartitions,1));

%cMeasures = { {'edit', '', ''}, {'','jaccard',''}, {'','rand',''}, {'','vi',''}, {'edgeRecon','',''}, {'pearsonRecon','',''}, {'', 'esVi',''},{'edgeRecon','vi','linear'}};
cMeasures = { {'', 'esVi',''} };
% measures by number of pairs of blockmodels
mDistances = zeros(size(cMeasures,2), size(cPartitions,1)-1);
    
% % construct sparse versions of the graphs in cSnapshots
% cSparseSnapshots = cell(size(cSnapshots,1), size(cSnapshots,2));
% for g = 1 : size(cSnapshots,2)
%    cSparseSnapshots{g} = sparse(cSnapshots{g}); 
% end

% number of vertices
vertNum = size(cSnapshots{1},1);

% construct vert-pos maps for the partitions
cVertPosMaps = cell(size(cPartitions,1), 1);
for s = 1 : size(cPartitions,1)
    % get number of partitions
    partNum = aPartInfo(s,3);
    currPart = 1;
    mVertPosMap = zeros(vertNum, partNum);
    for p = 1 : partNum
        vVertPos = cPartitions{s,p};
        for x = 1 : size(vVertPos,2)
            mVertPosMap(vVertPos(x), currPart) = 1.0; 
        end
        currPart = currPart + 1;
        cVertPosMaps{s} = mVertPosMap;
    end
end

% loop through each measure
for m = 1 : size(cMeasures,2)
    vMeasure = cMeasures{m}
    
    % loop through each consecutive graphs
    for t = 1 : (size(cSnapshots,2)-1)
%         t
% %         cSnapshots{t+1}
        vDist = compareOverlapBMNoFile(cVertPosMaps{t}, cSnapshots{t}, cVertPosMaps{t+1}, cSnapshots{t+1}, vMeasure{1}, vMeasure{2}, vMeasure{3},'temp.dis.csv',0.5);
        mDistances(m,t) = vDist;
    end
end % end of for of cMeasures



end