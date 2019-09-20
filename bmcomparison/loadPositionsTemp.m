function [mVertPosMap, aPartInfo] = loadPositionsTemp(sPosFilename, bHeader, bListFormat, vNum)
%
% Loads the set of positions.
%

aPartInfo = [];

% load the cPartitions file
fPosFilename = fopen(sPosFilename);
sLine = fgetl(fPosFilename);
currPart = 1;

while ischar(sLine)
    % we should be reading a header next
    if (bHeader)
        cInfo = textscan(sLine, '%d%d%d%f', 'Delimiter', ',');
        % we need to add one to time units
        startTime = cInfo{1} + 1;
        endTime = cInfo{2} + 1;
        partNum = cInfo{3};
        subseqCost = cInfo{4};
        
        aPartInfo = [startTime, endTime, partNum, subseqCost];
        
        % allocate size of mVertPosMap
        mVertPosMap = zeros(vNum, partNum);
        
        bHeader = false;
        
    else
        % not header, so partition information
        if (bListFormat)
            cVertices = textscan(sLine, '%d', 'Delimiter', ',');
            vVertPos = cVertices{1};
            for x = 1 : size(vVertPos,1)
               mVertPosMap(vVertPos(x), currPart) = 1.0; 
            end
            currPart = currPart + 1;
        else
            % read in partition
            cVertPos = textscan(sLine, '%f', 'Delimiter', ',');
            vVertPos = cVertPos{1};
            for x = 1 : floor((size(vVertPos,1)-1)/2)
                mVertPosMap(vVertPos(1), vVertPos(2*x)) = vVertPos(2*x + 1);
            end
        end
        
    end
    sLine = fgetl(fPosFilename);
end



end % end of function