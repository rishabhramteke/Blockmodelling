function [cSnapshots, cPartitions, aPartInfo] = loadMatpart(sBaseDir, sMatWildcard, sSubseqFile, bSparse, m, n, bAddOne)
%
% Loads the snapshots and subsequence information into cell arrays.
% For SeqiBloc, sequences of graphs and corresponding set of positions.
%
% bAddOne - whether to add 1 to vertex indices
%

% get the list of snapshot/matrices filenames
stMatFilenames = dir(fullfile(sBaseDir, sMatWildcard));

% initialise
cSnapshots = cell(1, size(stMatFilenames,1));

% load them one by ''one
for i = 1 : size(stMatFilenames,1)
    mSnapshot = [];
    % need to test if file is empty (Matlab doesn't check this in
    % csvread!!!
    fid = fopen(fullfile(sBaseDir, stMatFilenames(i).name));
    tline = fgetl(fid);
    commentLines = 0;
    bOkay = 1;
    if ~ischar(tline) || isempty(tline)    
        mSnapshot = zeros(m, n);
        % convert to sparse if all zeros
        if bSparse
            mSnapshot = sparse(mSnapshot);
        end
        bOkay = 0;
    elseif strcmp(tline(1), '#')
        % record the number of comment strings
        commentLines = commentLines + 1;
        tline = fgetl(fid);
        while ~isempty(tline) && strcmp(tline(1), '#')
            commentLines = commentLines + 1;
            tline = fgetl(fid);
        end
    end
    
    fclose(fid);
    
    addAmount = 0;
    if bAddOne
        addAmound = 1;
    end

    if bOkay   
        if bSparse
            mData = csvread(fullfile(sBaseDir, stMatFilenames(i).name), commentLines, 0);
            % no weights
            if size(mData,2) == 2
                mSnapshot = sparse(mData(:,1) +addAmound, mData(:,2) +addAmound, ones(size(mData,1),1), m, n);
            else
                % weighted
                mSnapshot = sparse(mData(:,1) + addAmound, mData(:,2) + addAmound, mData(:,3), m, n);
            end
        else
            mSnapshot = csvread(fullfile(sBaseDir, stMatFilenames(i).name), commentLines, 0);
        end
    end
    cSnapshots{i} = mSnapshot;
end % end of for

cPartitions = {};
aPartInfo = [];


% load the cPartitions file
fSubseqFile = fopen(sSubseqFile);
sLine = fgetl(fSubseqFile);
bHeader = true;
currPart = 1;
currSubseq = 0;
while ischar(sLine)
    % we should be reading a header next
    if bHeader
        cInfo = textscan(sLine, '%d%d%d%f', 'Delimiter', ',');
        % we need to add one to time units
        startTime = cInfo{1} + 1;
        endTime = cInfo{2} + 1;
        partNum = cInfo{3};
        subseqCost = cInfo{4};
        
        % create a new cPartition and cPartInfo entry
%         if size(cPartitions, 2) == 0
%             cPartitions = {cell(1,partNum)}
%         else
%             cPartitions = {cPartitions; cell(1,partNum)}
%         end
        
        vInfo = [startTime, endTime, partNum, subseqCost];
        if size(aPartInfo,2) == 0
            aPartInfo = vInfo;
        else
            aPartInfo = [aPartInfo ; vInfo];
        end
        
        bHeader = false;
        currSubseq = currSubseq + 1;
    else
        % not header, so partition information
              
        % read in partition
        cVertices = textscan(sLine, '%d', 'Delimiter', ',');
%         for y = 1 : size(cVertices{1})
%             cVertices{1}(y) = cVertices{1}(y) + 1;
%         end
        cPartitions{currSubseq,currPart} = cVertices{1}';

        % if true, this is the last partition to read before we expect
        % another header
        if currPart == partNum
            bHeader = true;
            currPart = 1;
        else
            currPart = currPart + 1;
        end
        
    end
    sLine = fgetl(fSubseqFile);
end



end % end of function


% function [mSnapshot] = loadEdgeList(sFilename, sDelimiter)
%     fSubseqFile = fopen(sSubseqFile);
%     sLine = fgetl(fSubseqFile);   
%     
%     bHeader = 1;
%     
%     while ischar(sLine)
%         % read header/comment lines
%         if bHeader
%             cLine = textscan(sLine, '%s', 'Delimiter', ' ');
%             if strcmp(cLine{1}, '#')
%                 sLine = fgetl(fSubseqFile);
%                 bHeader = 0
%                 continue;
%             end
%         end
%         
%         % no more header to read
%         cVertices = textscan(sLine, '%d', 'Delimiter', sDelimiter);
%             
%         sLine = fgetl(fSubseqFile);
%     end % end of while
% end % end of function