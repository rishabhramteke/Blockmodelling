function [mSparseGraph] = loadMatlabSparseGraph(sFilename, m, n, bMatlabFormat)
%
% Loads matlab sparse format graph.
%



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
    
fclose(fid);

if bOkay   
    if (bMatlabFormat)
        mData = csvread(sFilename, commentLines, 0);
        mSparseGraph = sparse(mData(:,1), mData(:,2), mData(:,3), m, n);
    else
        mData = csvread(sFilename, commentLines, 0);
        vWeights = ones(size(mData(:,1),1),1);
        mSparseGraph = sparse(mData(:,1), mData(:,2), vWeights, m, n);
    end
else
    mSparseGraph = sparse(m, n); 
end


 end % end of function