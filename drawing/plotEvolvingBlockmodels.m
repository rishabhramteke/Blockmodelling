function plotEvolvingBlockmodels(cmSnapshots, cPartitions, aPartInfo, subplotsPerRow)
%
% Plots a the snapshots specitied in cmSnapshots,
% as a sequence of adjacency matrices blockmodels, with the rows and columns
% rearranged according to cPartitions.  
%
% cmSnapshots - a cell vector of 2d matrices of the snapshots.
% cPartitions - a cell of vectors of the partitions (indices into the
% snapshots).  Each row is the partitions of one subsequence
% cMetaInfo - a cell of start index (into cmSnapshots), end index, 
% encoding costs
% aPartInfo - an array of partition info, each row is the information of one subsequence 
% subplotPerRow - number of subplots (snapshots) to draw per row
%

% make sure the sizes of the snapshots are the same.
assert(size(aPartInfo,1) == size(cPartitions,1));
% make sure the length of the subsequences is the same as the number of
% snapshots
assert(size(cmSnapshots,2) == sum(aPartInfo(end,2)));

% plot
figure;

% calculate the size of the subplot grid
totalRow = ceil(size(cmSnapshots,2) / subplotsPerRow);
totalCol = subplotsPerRow;

currRow = 0;
currCol = 1;


for currSubseq = 1 : size(aPartInfo,1)
    % get subsequence information
    partNum = aPartInfo(currSubseq,3);
    startTime = aPartInfo(currSubseq,1);
    endTime = aPartInfo(currSubseq,2);
    
    % min values of each partition
    vPartMinVal = zeros(1,partNum);
    
    % loop throuch each partition of the subsequence
    for currPart = 1 : partNum
        % sort each partition and store its minimum values
        cPartitions{currSubseq, currPart} = sort(cPartitions{currSubseq, currPart});
        vPartMinVal(currPart) = min(cPartitions{currSubseq, currPart});
    end % end of for looping partitions        
    
    % sort the order
    [vSortedOrder, vOrderIndex] = sort(vPartMinVal); 
    
    if currSubseq == 1
        % sort the order of the partitions based on their minimum value
        cPartitions(currSubseq, [1:1:currPart]) = cPartitions(currSubseq, vOrderIndex);
    else
        % sort the ordering of the partitoins based on the previous
        % ordering
        vOrderIndexBasedOnPrev = matchOrdering(cPartitions(currSubseq-1, [1:1:currPart]), cPartitions(currSubseq, [1:1:currPart]));
        cPartitions(currSubseq, [1:1:currPart]) = cPartitions(currSubseq, vOrderIndexBasedOnPrev);
    end
    % boundary locations
    vBoundary = zeros(1, partNum);
    vI = [];
    startLocation = 1;
    for j = 1 : partNum% sort the order of the partitions based on their minimum value
        vI = [vI cPartitions{currSubseq, j}];
        vBoundary(j) = startLocation + size(cPartitions{currSubseq, j},2);
        startLocation = startLocation + size(cPartitions{currSubseq,j},2);
    end
    
    
    for t = 1 : endTime - startTime + 1
        % subplot variables
        if currCol > subplotsPerRow
            currRow = currRow + 1;
            currCol = 1;
        end        
    
        subplot(totalRow, totalCol, currRow * subplotsPerRow + currCol);
        % get the appropriate snapshot
        mAdj = cmSnapshots{currRow * subplotsPerRow + currCol};
        % plot rearranged snapshot
        drawBlockmodel(mAdj(vI, vI), vBoundary, vBoundary); 
        %drawBlockmodel(mAdj(vI, vI), [], []); 
%         drawBlockmodel(mAdj, [], []);
    
        currCol = currCol + 1;            
    end    
    
end % end of for


end % end of function


function [vNewOrder] = matchOrdering(cP1Order, cP2)
%
% Computes the indices that the partitions in cP2 should be reordred to.
%

mJacVals = zeros(size(cP1Order,2), size(cP2,2));
for p1 = 1 : size(cP1Order,2)
   for p2 = 1 : size(cP2,2)
       mJacVals(p1,p2) = calcJaccard(cP1Order{p1}, cP2{p2});
   end
end

vNewOrder = [];
for p1 = 1 : size(cP1Order,2)
    % find all jaccard values that are greater than zero for p1
%     vNonZeroIndices = find(mJacVals(p1,:) > 0);
%     if size(vNonZeroIndices,1) == 0
%         % random assign an unallocated one
%         for j = 1 : size(
%         vNewOrder = [vNewOrder 
    
    vIndices = mJacVals(p1,:);
    % sort them by jaccard values
%     [vSorted, vSortedIndices] = sort(mJacVals(p1,vNonZeroIndices), 'descend');
    [vSorted, vSortedIndices] = sort(vIndices, 'descend');
    
    % reordered non zero indices according to sortvSorted
%     vReordered = vNonZeroIndices(vSortedIndices);
    vReordered = vSortedIndices;
    for i = 1 : size(vReordered, 2)
        % only add to vNewOrder if it hasn't been added yet
        if size(find(vNewOrder == vReordered(i)),2) == 0
            vNewOrder = [vNewOrder vReordered(i)];
            break;
        end
        
    end
end



end % end of function



function [jaccardVal] = calcJaccard(vP1, vP2)
%
% Compute the jaccard of the two sets vP1 and vP2.
%

% vP1
% vP2
jaccardVal = size(intersect(vP1, vP2),2) / size(union(vP1, vP2),2);
end % end of function