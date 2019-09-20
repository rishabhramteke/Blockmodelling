function [AR,RI,MI,HI] = edgeSetRandIndex(vMembership1, mAdj1, vMembership2, mAdj2)
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
% vEdgeBlock1 = zeros(unionEdgeNum,1);
% vEdgeBlock2 = zeros(unionEdgeNum,1);
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
% % call randIndex
% [AR,RI,MI,HI] = randIndex(vEdgeBlock1, vEdgeBlock2);


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
% mInterAdj = mAdj1 & mAdj2;

% number of positions
posNum1 = size(unique(vMembership1), dim);
posNum2 = size(unique(vMembership2), dim);


mPosOverlap = zeros(posNum1, posNum2);
for v1 = 1 : size(vMembership1,dim)
    mPosOverlap(vMembership1(v1), vMembership2(v1)) = mPosOverlap(vMembership1(v1), vMembership2(v1)) + 1;
end


blockNum1 = posNum1 * posNum1;
blockNum2 = posNum2 * posNum2;

mInter = zeros(blockNum1, blockNum2);

% find all edges that appear in one graph but not the other
[vRow, vCol] = find(mXorAdj);
for e = 1 : size(vRow,1)
    mInter(vMembership1(vRow(e)), vMembership1(vCol(e))) =...
        mInter(vMembership1(vRow(e)), vMembership1(vCol(e))) + 1;
end



% now we convert mInterDistrib to the number of common edges and non-edges
for r = 1 : size(mInterDistrib, 1)
    for c = 1 : size(mInterDistrib, 2)
        posRowIndex1 = floor(r / (posNum1)) + 1;
        posColIndex1 = floor(mod(r,posNum1));
        if (posColIndex1 == 0)
            posRowIndex1 = posRowIndex1 - 1;
            posColIndex1 = posColIndex1 + posNum1;
        end
        posRowIndex2 = floor(c / (posNum2)) + 1;
        posColIndex2 = floor(mod(c,posNum2));
        if (posColIndex2 == 0)
            posRowIndex2 = posRowIndex2 - 1;
            posColIndex2 = posColIndex2 + posNum2;
        end        
        if (mPosOverlap(posRowIndex1,posRowIndex2) > 0 && mPosOverlap(posColIndex1,posColIndex2) > 0)
            mInter(r,c) =...
                (mPosOverlap(posRowIndex1,posRowIndex2) * mPosOverlap(posColIndex1,posColIndex2) - mInter(r,c));
        end
    end
end






n=sum(sum(mInter));
nis=sum(sum(mInter,2).^2);		%sum of squares of sums of rows
njs=sum(sum(mInter,1).^2);		%sum of squares of sums of columns

t1=nchoosek(n,2);		%total number of pairs of entities
t2=sum(sum(mInter.^2));	%sum over rows & columnns of nij^2
t3=.5*(nis+njs);

%Expected index (for adjustment)
nc=(n*(n^2+1)-(n+1)*nis-(n+1)*njs+2*(nis*njs)/n)/(2*(n-1));

A=t1+t2-t3;		%no. agreements
D=  -t2+t3;		%no. disagreements

if t1==nc
   AR=0;			%avoid division by zero; if k=1, define Rand = 0
else
   AR=(A-nc)/(t1-nc);		%adjusted Rand - Hubert & Arabie 1985
end

RI=A/t1;			%Rand 1971		%Probability of agreement
MI=D/t1;			%Mirkin 1970	%p(disagreement)
HI=(A-D)/t1;	%Hubert 1977	%p(agree)-p(disagree)





end

function [blockNum] = posPair2BlockNum(rowPos, colPos, posNum)

blockNum = (rowPos - 1) * posNum + colPos;

end