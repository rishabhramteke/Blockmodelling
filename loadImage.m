function [mImage] = loadImage(sImageFilename, posNum)
%
% Loads the set of positions.
%


% load the cPartitions file
fImageFilename = fopen(sImageFilename);
sLine = fgetl(fImageFilename);
currLine = 1;
mImage = zeros(posNum, posNum);

while ischar(sLine)
    cLine = textscan(sLine, '%f', 'Delimiter', ',');
    vEntries = cLine{1};
    
    mImage(currLine,:) = vEntries;
    
    currLine = currLine + 1;
        
    sLine = fgetl(fImageFilename);
end



end % end of function