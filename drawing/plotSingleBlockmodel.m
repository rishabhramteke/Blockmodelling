function plotSingleBlockmodel(mSnapshot, mVertPosMap, sTitle, bUseVertIndices, bSortPos)
%
% Plots the snapshot specified in mSnapshots,
% with the rows and columns rearranged according to mVertPosMap.  
%
% mSnapshot - adjacency matrix.
% mVertPosMap - membershiop matrix.
% sTitle - title of figure.
% bUseVertIndices   - Use the vertex indices as markings on image.
%
% @author Jeffrey Chan, 2012
%


% plot
figure;
title(sTitle);

partNum = size(mVertPosMap,2);
cPartitions = cell(partNum,1);
% intialialise cPartiitons
for p = 1 : partNum
    cPartitions{p} = [];
end

% convert to a cell of partitions
for v = 1 : size(mVertPosMap,1)
    part = find(mVertPosMap(v,:));
    cPartitions{part} = [cPartitions{part} v];
end

% if we are missing a part, we delete it
vNoGood = [];
for p = 1 : partNum
    if size(cPartitions{p}) == 0
        vNoGood = [p vNoGood];
    end
end

for i = 1 : length(vNoGood)
    cPartitions(vNoGood(i)) = [];
    partNum = partNum - 1;
end
    
% min values of each partition
vPartMinVal = zeros(1,partNum);
    
% loop throuch each partition of the subsequence
for currPart = 1 : partNum
    % sort each partition and store its minimum values
    cPartitions{currPart} = sort(cPartitions{currPart});
    vPartMinVal(currPart) = min(cPartitions{currPart});
end % end of for looping partitions        

% sort the order
if bSortPos
    [vSortedOrder, vOrderIndex] = sort(vPartMinVal); 

    % sort the order of the partitions based on their minimum value
    cPartitions = cPartitions(vOrderIndex);
end



% boundary locations

    vBoundary = zeros(1, partNum-1);
    vI = [];
    startLocation = 1;
    for j = 1 : partNum
        vI = [vI cPartitions{j}];
    end
    for j = 1 : partNum-1
        vBoundary(j) = startLocation + size(cPartitions{j},2)-0.5;
        startLocation = startLocation + size(cPartitions{j},2); 
    end



    % plot rearranged snapshot
    if bUseVertIndices
        drawBlockmodel(mSnapshot(vI, vI), vBoundary, vBoundary, vI, vI); 
    else
        drawBlockmodel(mSnapshot(vI, vI), vBoundary, vBoundary); 
    end
    


end % end of function
