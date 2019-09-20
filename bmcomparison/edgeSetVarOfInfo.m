function [distance] = edgeSetVarOfInfo(vMembership1, mAdj1, vMembership2, mAdj2, bEdit)
%
% Computes the variation of information for edge sets overlap over the blocks 
% of the two blockmodels (represented by vMembership vectors and the adjacency
% matrices they were derived from).
%
% vMembership1 - vector of membership (each element is 1..K, where K is the
% number of clusters in clustering set 1)
% mAdj1 - Adjacency matrix 1, in sparse format.
% vMembership2 - vector of membership (each element is 1..K', where K' is the
% number of clusters in clustering set 2)
% mAdj1 - Adjacency matrix 2, in sparse format.
% bEdit - Whether an edit distance.
%



% [r,c] = size(vMembership1);
% dim = 1;
% if c > r
%     dim = 2;
% end
% 
% % make sure the adjacency matrices are sparse
% assert(issparse(mAdj1) && issparse(mAdj2));
% assert(size(vMembership1,dim) == size(vMembership2,dim));
% 
% % map each edge in the union of the two matrices to an integer 
% %(continuous integer range)
% % the edge mapping to an index doesn't have to be the same for the two
% % matrices, as we only care about the distributions, not the actual
% % matchings in VI
% unionEdgeNum = size(find(mAdj1 | mAdj2),1);
% mInterAdj = mAdj1 & mAdj2;
% mIn1Not2 = mAdj1 - mInterAdj;
% mIn2Not1 = mAdj2 - mInterAdj;
% 
% 
% 
% % number of positions
% vPosNum1 = size(unique(vMembership1), dim);
% vPosNum2 = size(unique(vMembership2), dim);
% 
% 
% % construct membership vectors for this mapping
% vEdgeBlock1 = zeros(1,unionEdgeNum);
% vEdgeBlock2 = zeros(1,unionEdgeNum);
% 
% % find all edges in both adjacency matrices
% [vRow, vCol] = find(mInterAdj);
% for e = 1 : size(vRow,1)
%     % look up positions for each edge
%     blockNum1 = posPair2BlockNum(vMembership1(vRow(e)), vMembership1(vCol(e)), vPosNum1);
%     blockNum2 = posPair2BlockNum(vMembership2(vRow(e)), vMembership2(vCol(e)), vPosNum2);
%     vEdgeBlock1(e) = blockNum1;
%     vEdgeBlock2(e) = blockNum2;
% end
% 
% % now add the positions that are in 1 but not in 2 as singletons
% currBlockNum2 = vPosNum2 * vPosNum2 + 1;
% [vIn1Not2Row, vIn1Not2Col] = find(mIn1Not2);
% for ee = 1 : size(vIn1Not2Row)
%    blockNum1 = posPair2BlockNum(vMembership1(vIn1Not2Row(ee)), vMembership1(vIn1Not2Col(ee)), vPosNum1);
%    vEdgeBlock1(ee + size(vRow,1)) = blockNum1;
%    vEdgeBlock2(ee + size(vRow,1)) = currBlockNum2;
%    currBlockNum2 = currBlockNum2 + 1;
% end
% 
% % now add the positions that are in 2 but not in 1 as singletons
% currBlockNum1 = vPosNum1 * vPosNum1 + 1;
% [vIn2Not1Row, vIn2Not1Col] = find(mIn2Not1);
% for ee = 1 : size(vIn2Not1Row)
%    blockNum2 = posPair2BlockNum(vMembership2(vIn2Not1Row(ee)), vMembership2(vIn2Not1Col(ee)), vPosNum1);
%    vEdgeBlock2(ee + size(vRow,1) + size(vIn1Not2Row,1)) = blockNum2;
%    vEdgeBlock1(ee + size(vRow,1) + size(vIn1Not2Row,1)) = currBlockNum1;
%    currBlockNum1 = currBlockNum1 + 1;
% end
% 
% 
% distance = varOfInfo(vEdgeBlock1, vEdgeBlock2);


[r,c] = size(vMembership1);
dim = 1;
if c > r
    dim = 2;
end

% make sure the adjacency matrices are sparse
assert(issparse(mAdj1) && issparse(mAdj2));
assert(size(vMembership1,dim) == size(vMembership2,dim));

% find all edges in one matrix but not in the other matrix
mXorAdj = xor(mAdj1,mAdj2);
mInterAdj = mAdj1 & mAdj2;
mIn1Not2 = mAdj1 - mInterAdj;
mIn2Not1 = mAdj2 - mInterAdj;

% number of positions
posNum1 = size(unique(vMembership1), dim);
posNum2 = size(unique(vMembership2), dim);
vPosSize1 = zeros(1, posNum1);
vPosSize2 = zeros(1, posNum2);
elemNum = size(vMembership1,dim) * size(vMembership1, dim);

% find position sizes
for p = 1 : posNum1
    vPosSize1(p) = size(find(vMembership1 == p), dim);
end
for p = 1 : posNum2
    vPosSize2(p) = size(find(vMembership2 == p), dim);
end


mPosOverlap = sparse(posNum1, posNum2);
for v1 = 1 : size(vMembership1,dim)
    mPosOverlap(vMembership1(v1), vMembership2(v1)) = mPosOverlap(vMembership1(v1), vMembership2(v1)) + 1;
end


blockNum1 = posNum1 * posNum1;
blockNum2 = posNum2 * posNum2;
% vBlockSize1 = sparse(1, blockNum1);
% vBlockSize1 = sparse(1, elemNum);
% vBlockSize2 = sparse(1, elemNum);
% % find block sizes
% for pr = 1 : posNum1
%     for pc = 1 : posNum1
%         vBlockSize1((pr-1) * posNum1 + pc) = vPosSize1(pr) * vPosSize1(pc);
%     end
% end
% 
% for pr = 1 : posNum2
%     for pc = 1 : posNum2
%         vBlockSize2((pr-1) * posNum2 + pc) = vPosSize2(pr) * vPosSize2(pc);
%     end
% end

mBlockSize1 = sparse(elemNum, elemNum);
mBlockSize2 = sparse(elemNum, elemNum);
% find block sizes
for pr = 1 : posNum1
    for pc = 1 : posNum1
        mBlockSize1(pr, pc) = vPosSize1(pr) * vPosSize1(pc);
    end
end

for pr = 1 : posNum2
    for pc = 1 : posNum2
        mBlockSize2(pr, pc) = vPosSize2(pr) * vPosSize2(pc);
    end
end




% mInterDistrib = sparse(blockNum1, blockNum2);
mInterDistrib = sparse(elemNum, elemNum);

% for all block-block that have a non-zero overlap, we initialise with the
% maximum amount of overlap
[vRow, vCol, vVal] = find(mPosOverlap);
for p1 = 1 : size(vRow,1)
    for p2 = 1 : size(vRow,1)
  
        % convert to mInterDistrib indexing
        rowOverlapIndex = (vRow(p1) - 1) * posNum1 + vRow(p2);
        colOverlapIndex = (vCol(p1) - 1) * posNum2 + vCol(p2);
        mInterDistrib(rowOverlapIndex, colOverlapIndex) = mPosOverlap(vRow(p1), vCol(p1)) * mPosOverlap(vRow(p2), vCol(p2));
    end
end


% NEW addition
% NON-EDIT VERSION
if ~bEdit 
    % match mIn1Not2 and mIn2Not1 across blocks.  Crucial part is we add 1
    % for mIn1Not2 and subtract 1 for mIn2Not1
    mBlockMatchings = zeros(blockNum1, blockNum2);
   
    [vRow1not2, vCol1not2] = find(mIn1Not2);
    for e = 1 : size(vRow1not2,1)
        rowIndex = (vMembership1(vRow1not2(e)) - 1) * posNum1 + vMembership1(vCol1not2(e));
        colIndex = (vMembership2(vRow1not2(e)) - 1) * posNum2 + vMembership2(vCol1not2(e));
        mBlockMatchings(rowIndex, colIndex) = mBlockMatchings(rowIndex, colIndex) + 1;
    end
    
    [vRow2not1, vCol2not1] = find(mIn2Not1);
    for e = 1 : size(vRow2not1,1)
        rowIndex = (vMembership1(vRow2not1(e)) - 1) * posNum1 + vMembership1(vCol2not1(e));
        colIndex = (vMembership2(vRow2not1(e)) - 1) * posNum2 + vMembership2(vCol2not1(e));        
        mBlockMatchings(rowIndex, colIndex) = mBlockMatchings(rowIndex, colIndex) - 1;
    end  
    
    
    % find non-zero entries and insert appropriate singleton entries
    currCol = blockNum2 + 1;
    currRow = blockNum1 + 1;
    [vRowMatch, vColMatch, vMatchVal] = find(mBlockMatchings);
    for b = 1 : size(vRowMatch,1) 
        while vMatchVal(b) > 0
            mInterDistrib(vRowMatch(b), vColMatch(b)) = mInterDistrib(vRowMatch(b), vColMatch(b)) - 1;
            % add new singleton entry (column) for each different one
            mInterDistrib(vRowMatch(b), currCol) = 1;
            % update blocksize
            [rowPos, colPos] = convertIndex(vColMatch(b), posNum2);
            [newRowPos, newColPos] = convertIndex(currCol, posNum2);
            mBlockSize2(rowPos, colPos) = mBlockSize2(rowPos, colPos) - 1;
            mBlockSize2(newRowPos, newColPos) = 1;    
            currCol = currCol + 1;            
            vMatchVal(b) = vMatchVal(b) - 1;
        end
        while vMatchVal(b) < 0
            mInterDistrib(vRowMatch(b), vColMatch(b)) = mInterDistrib(vRowMatch(b), vColMatch(b)) - 1;
            % add new singleton entry (row) for each different one
            mInterDistrib(currRow, vColMatch(b)) = 1;
            % update blocksize
            [rowPos, colPos] = convertIndex(vRowMatch(b), posNum1);
            [newRowPos, newColPos] = convertIndex(currRow, posNum1);
            mBlockSize1(rowPos, colPos) = mBlockSize1(rowPos, colPos) - 1;
            mBlockSize1(newRowPos, newColPos) = 1;     
            currRow = currRow + 1;
            vMatchVal(b) = vMatchVal(b) + 1;
        end
    end
else
    % EDIT DISTANCE VERSION
    % find all edges that appear in one graph but not the other
    % [vRow, vCol] = find(mXorAdj);
    currCol = blockNum2 + 1;
    [vRow1not2, vCol1not2] = find(mIn1Not2);
    for e = 1 : size(vRow1not2,1)
        rowIndex = (vMembership1(vRow1not2(e)) - 1) * posNum1 + vMembership1(vCol1not2(e));
        colIndex = (vMembership2(vRow1not2(e)) - 1) * posNum2 + vMembership2(vCol1not2(e));
        mInterDistrib(rowIndex, colIndex) = mInterDistrib(rowIndex, colIndex) - 1;
        % add new singleton entry (column) for each different one
        mInterDistrib(rowIndex, currCol) = 1;
        % update blocksize
        [rowPos, colPos] = convertIndex(colIndex, posNum2);
        [newRowPos, newColPos] = convertIndex(currCol, posNum2);
        mBlockSize2(rowPos, colPos) = mBlockSize2(rowPos, colPos) - 1;
        mBlockSize2(newRowPos, newColPos) = 1;    
        currCol = currCol + 1;
    end


    [vRow2not1, vCol2not1] = find(mIn2Not1);
    currRow = blockNum1 + 1;
    for e = 1 : size(vRow2not1,1)
        rowIndex = (vMembership1(vRow2not1(e)) - 1) * posNum1 + vMembership1(vCol2not1(e));
        colIndex = (vMembership2(vRow2not1(e)) - 1) * posNum2 + vMembership2(vCol2not1(e));
        mInterDistrib(rowIndex, colIndex) = mInterDistrib(rowIndex, colIndex) - 1;
        % add new singleton entry (row) for each different one
        mInterDistrib(currRow, colIndex) = 1;
        % update blocksize
        [rowPos, colPos] = convertIndex(rowIndex, posNum1);
        [newRowPos, newColPos] = convertIndex(currRow, posNum1);
        mBlockSize1(rowPos, colPos) = mBlockSize1(rowPos, colPos) - 1;
        mBlockSize1(newRowPos, newColPos) = 1;     
        currRow = currRow + 1;
    end
end








% normalise
mInterDistrib = mInterDistrib ./ elemNum;


% now we convert mInterDistrib to the number of common edges and non-edges
% for r = 1 : size(mInterDistrib, 1)
%     for c = 1 : size(mInterDistrib, 2)
%         posRowIndex1 = floor(r / (posNum1)) + 1;
%         posColIndex1 = floor(mod(r,posNum1));
%         if (posColIndex1 == 0)
%             posRowIndex1 = posRowIndex1 - 1;
%             posColIndex1 = posColIndex1 + posNum1;
%         end
%         posRowIndex2 = floor(c / (posNum2)) + 1;
%         posColIndex2 = floor(mod(c,posNum2));
%         if (posColIndex2 == 0)
%             posRowIndex2 = posRowIndex2 - 1;
%             posColIndex2 = posColIndex2 + posNum2;
%         end        
%         if (mPosOverlap(posRowIndex1,posRowIndex2) > 0 && mPosOverlap(posColIndex1,posColIndex2) > 0)
%             mInterDistrib(r,c) =...
%                 (mPosOverlap(posRowIndex1,posRowIndex2) * mPosOverlap(posColIndex1,posColIndex2) - mInterDistrib(r,c)) / elemNum;
%         end
%     end
% end


mDistrib1 = mBlockSize1 ./ elemNum;
mDistrib2 = mBlockSize2 ./ elemNum;


mutualVal = 0;
[vInterRow, vInterCol, vInterVals] = find(mInterDistrib);

for m = 1 : size(vInterRow)
% for i = 1 : size(mInterDistrib, 1)
%     for j = 1 : size(mInterDistrib, 2)
        % mylog returns negative values, so we adding negative values, so
        % we need to negate it first
    r = vInterRow(m);
    c = vInterCol(m);
    [rowRPos, rowCPos] = convertIndex(r, posNum1);
    [colRPos, colCPos] = convertIndex(c, posNum2);    
    if (mDistrib1(rowRPos, rowCPos) * mDistrib2(colRPos,colCPos) == 0)
        rowRPos
        rowCPos
        colRPos
        colCPos
    end
    assert(mDistrib1(rowRPos, rowCPos) * mDistrib2(colRPos,colCPos) > 0);
    mutualVal = mutualVal + mInterDistrib(r,c) * myLog(mInterDistrib(r,c) / (mDistrib1(rowRPos, rowCPos) * mDistrib2(colRPos,colCPos)));
%     end
end

entropyVal1 = myEntropyMat(mDistrib1);
entropyVal2 = myEntropyMat(mDistrib2);

distance = entropyVal1 + entropyVal2 - 2 * mutualVal;

end



function [value] = myEntropyMat(mDistribution)
%
% compute entropy
%


value = 0;
[vR, vC, vV] = find(mDistribution);
for p = 1 : size(vR)
    value = value + mDistribution(vR(p), vC(p)) * log2(mDistribution(vR(p), vC(p)));
end

value = -value;

end


function [posRowIndex, posColIndex] = convertIndex(index, posNum)
%
% convert mInterDistrib indexing to position indexing
%

posRowIndex = floor(index / (posNum)) + 1;
posColIndex = floor(mod(index,posNum));
if (posColIndex == 0)
    posRowIndex = posRowIndex - 1;
    posColIndex = posColIndex + posNum;
end


end


% function [blockNum] = posPair2BlockNum(rowPos, colPos, posNum)
% 
% blockNum = (rowPos - 1) * posNum + colPos;
% 
% end