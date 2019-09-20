%Compute the cost of individual snapshot encoding
function [c, c3] = encodingCost(adjMatrix, mMembership)
    
    vertNum = size(mMembership,1);
    % get the indices for each block
    posNum = size(mMembership,2);
    cvPosIndices = cell(1, posNum);
    
    mNewMembership = discretise(mMembership);
    for p = 1 : posNum
        cvPosIndices{p} = find(mNewMembership(:,p));
    end    

    c1 = vertNumCost(vertNum);
    c2 = posCost(cvPosIndices, vertNum);
    c3 = snapshotCost(adjMatrix, cvPosIndices);
 
    c = c1 + c2 + c3;
%     fprintf('%.3f, %.3f, %.3f, %.3f\n', c1, c2, c3, c);
end

%Compute the cost to encode the number of vertices
function c1 = vertNumCost(n)
    c1 = log2(log2(n));
end

%Compute the cost to encode positions
function c2 = posCost(cvPosIndices, vertNum)
    posNum = length(cvPosIndices);

    entropyVal = 0;
    for p = 1 : posNum
        prob = length(cvPosIndices{p}) / vertNum;
        entropyVal = entropyVal - prob * mylog2(prob);
    end
    
    if posNum > 1
        c2 = log2(log2(posNum)) + vertNum * entropyVal;
    else
        c2 = 1 + vertNum * entropyVal;
    end
end



%Compute the cost to encode the snapshot
function c3 = snapshotCost(mAdj, cvPosIndices)

    posNum = length(cvPosIndices);
    c3 = 0;
    
    for r = 1 : posNum
        for c = 1 : posNum
            mSubAdj = mAdj(cvPosIndices{r}, cvPosIndices{c});
            subAdjLen = size(mSubAdj,1) * size(mSubAdj,2);
            if subAdjLen ~= 0
                oneCount = nnz(mSubAdj);
                oneProb = oneCount / subAdjLen;
                zeroProb = 1 - oneProb;
                entropyVal = -zeroProb * mylog2(zeroProb) - oneProb*mylog2(oneProb);
                if oneCount > 1
                    c3  = c3 + log2(log2(oneCount)) + subAdjLen * entropyVal;
                else
                    c3  = c3 + 1 + subAdjLen * entropyVal;
                end
            else
                c3 = c3 + 1;
            end
        end
    end
    

end 

function [val] = mylog2(prob)
    if prob == 0
        val = 0;
    else
        val = log2(prob);
    end
end

