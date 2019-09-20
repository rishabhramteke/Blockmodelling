function [mSparseGraph, mDim, nDim] = loadMatlabGraph(sFilename, m, n, bMatlabFormat, bSparse, bAddOneToIndex, bHeader)
%
% Loads matlab format graph.
%
% bSparse -         Should be true if input graph is in sparse format
%                   (i.e., edgelist).
% bMatlabFormat -   Edge list in (u, v, w) format.
% bAddOneToIndex -  Add one to the indexing of the vertices (input file
%                   starts indexing at 0).
% bHeader -         Whether to ignore first line as the "header" line.
%

mDim = m;
nDim = n;


% need to test if file is empty (Matlab doesn't check this in csvread!!!
fid = fopen(sFilename);
tline = fgetl(fid);
commentLines = 0;
bOkay = 1;
if ~ischar(tline) || isempty(tline)
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

% see if we need to read a header line
if bOkay && ~isempty(tline) && bHeader
    commentLines = commentLines + 1;
    % read in the size of the graph
    stSize = regexp(tline, '\s+','split');
    assert(length(stSize) == 2);
    mDim = str2double(stSize{1});
    nDim = str2double(stSize{2});
    % no need to read in next line
end

    
fclose(fid);

% everything okay (hopefully)
if bOkay   
    if (bSparse)
        if (bMatlabFormat)
            mData = csvread(sFilename, commentLines, 0);
            if (bAddOneToIndex)
                mSparseGraph = sparse(mData(:,1)+1, mData(:,2)+1, mData(:,3), mDim, nDim);
            else
                mSparseGraph = sparse(mData(:,1), mData(:,2), mData(:,3), mDim, nDim);
            end
        else
            mData = csvread(sFilename, commentLines, 0);
            vWeights = ones(size(mData(:,1),1),1);
            if (bAddOneToIndex)
                mSparseGraph = sparse(mData(:,1)+1, mData(:,2)+1, vWeights, mDim, nDim);
            else
                mSparseGraph = sparse(mData(:,1), mData(:,2), vWeights, mDim, nDim);
            end
        end
    else
        % not sparse 
        mSparseGraph = sparse(csvread(sFilename, commentLines, 0));
    end
else
    % empty graph, we just construct a empty graph
    mSparseGraph = sparse(m, n); 
end


 end % end of function