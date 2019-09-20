function [denDis] = bmOverlapDenDist(mMembership1, mData1, mMembership2, mData2,...
    sReconDist)
%
% Computes the density distances between clustering/bm 1 and
% 2.
%
% sReconDist - Reconstruction distance to use.
%

% construct image matrix (assume partitional for now)
mImage1 = constructImage(mMembership1, mData1);
mImage2 = constructImage(mMembership2, mData2);

switch sReconDist
    case 'edgeRecon'
        % compute density distance
        denDis = bmCompareRecon(mMembership1, mMembership2, mImage1, mImage2);
    case 'edgeReconEqual'
        % compute weight adjusted measure.  (not, mData1 must equal mData2)
        assert(nnz(mData1 == mData2) == size(mData1,1) * size(mData2,2));
        expectedDensity = (sum(sum(mData1)) / size(mData1,1)^2);
        mAdjWeight = mData1 - expectedDensity;
        denDis = bmCompareRecon(mMembership1, mMembership2, mImage1, mImage2, mAdjWeight);        
    case 'edgeReconDeg'
        % compute weight adjusted measure.  (not, mData1 must equal mData2)
        assert(nnz(mData1 == mData2) == size(mData1,1) * size(mData2,2));
        vInDeg = sum(mData1,1);
        vOutDeg = sum(mData1,2);
        % outer product
        edgeNum = sum(sum(mData1));
        mExpectedDensity = vInDeg * vOutDeg / (edgeNum^2);
        mAdjWeight = mData1 - mExpectedDensity;       
        denDis = bmCompareRecon(mMembership1, mMembership2, mImage1, mImage2, mAdjWeight);
    case 'blockRecon'
        % block distance
        
    case 'pearsonRecon'
        denDis = bmPearsonLabel(mMembership1, mMembership2, mImage1, mImage2, mData1, mData2);
    case 'roleRecon'
        denDis = bmRoleLabel(mMembership1, mMembership2, mImage1, mImage2);
    case 'edit'
        % compute edit distance
        mXor = xor(mData1, mData2);
        denDis = size(find(mXor),1);
    case 'edgeReconKL'
        denDis = bmCompareReconKL(mMembership1, mMembership2, mImage1, mImage2);
    case 'edgeReconHell'
        denDis = bmCompareReconFDiv(mMembership1, mMembership2, mImage1, mImage2, 'hellinger');
    case 'edgeReconPearson'
        denDis = bmCompareReconFDiv(mMembership1, mMembership2, mImage1, mImage2, 'pearson');
    otherwise
        warning('Unknown sReconDist = %s specified', sReconDist);
        return;
end


end % end of function




