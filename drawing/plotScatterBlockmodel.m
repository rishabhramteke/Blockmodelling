function [vXcopy, vYcopy] = plotScatterBlockmodel(cSnapshots, cPartitions, partNum, vPartOrder)
%
% Plots the union graph of the snapshots specified in cmSnapshots,
% with the rows and columns rearranged according to cPartitions.  
%
% Calls scatterplot to draw the actual blockmodels.
%
% Intially used for drawing larger blockmodels like the BGP ones.
%
% cmSnapshots - a cell vector of 2d matrices of the snapshots.
% cPartitions - a cell of vectors of the partitions (indices into the
% snapshots).  Each row is the partitions of one subsequence
% aPartInfo - an array of partition info, each row is the information of one subsequence 
% subplotPerRow - number of subplots (snapshots) to draw per row
%


vX = [];
vY = [];
% compute the union of cSnapshots
for s = 1 : size(cSnapshots,2)
    currAdj = cSnapshots{s};
    if s == 1
        vX = currAdj(:,1);
        vY = currAdj(:,2);
    else
        for x = 1 : size(currAdj(:,1))
            if isempty(any(vX == currAdj(x,1))) && isempty(any(vY == currAdj(x,2))),
                    vX = [vX; currAdj(x,1)];
                    vY = [vY; currAdj(x,2)];
            end
        end
    end
end

% boundary locations
vBoundary = zeros(1, partNum);
vI = [];
startLocation = 1;
for j = 1 : partNum% sort the order of the partitions based on their minimum value
    vI = [vI cPartitions{7, vPartOrder(j)}];
    vBoundary(j) = startLocation + size(cPartitions{7, vPartOrder(j)},2);
    startLocation = startLocation + size(cPartitions{7,vPartOrder(j)},2);
%     vI = [vI cPartitions{1, j}];
%     vBoundary(j) = startLocation + size(cPartitions{1, j},2);
%     startLocation = startLocation + size(cPartitions{1,j},2);
end

vXcopy = vX;
vYcopy = vY;

for i = 1 : size(vXcopy)
    vNewLoc = find(vI == vXcopy(i));
    if isempty(vNewLoc)
        display "empty"
    end
end


for i = 1 : size(vYcopy)
    vNewLoc = find(vI == vYcopy(i));
    if isempty(vNewLoc)
        display "empty"
    end
    %vYcopy(i) = vNewLoc;
    if vNewLoc > 14000
        vYcopy(i) = vNewLoc - 3600;
    else
        vYcopy(i) = vNewLoc;
    end    
end

vvXcopy = [vXcopy ; vYcopy];
vvYcopy = [vYcopy ; vXcopy];

% plot rearranged snapshot
scatterBlockmodel(vvXcopy,vvYcopy, vBoundary, vBoundary); 




end % end of function