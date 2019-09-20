function [dist] = bmRoleLabel(mMembership1, mMembership2, mImage1, mImage2, mData, sNullModel, gamma)
% 
% Compute the reconstuction distance based on the Reichardt objective.
% We use the error formulation rather than the quality function.
% We assume we are working on the same adjacency matrix (TODO: to see if
% the link/non-link matching weights can work when adjacency matrices are
% different.
%

% determine how to compute mP, the null model probabilities
mP = [];

switch sNullModel
    case 'equal'
        totalEdgeWeight = sum(sum(mData));
        mP = ones(size(mData,1)) * (totalEdgeWeight / size(mData,1)^2);
    case 'power'
        % out degree
        vKout = sum(mData,2);
        vKin = sum(mData,1);
        edgeNum = sum(sum(find(mData > 0));
        mP = vKout * vKin' / edgeNum;
    otherwise
        warning('Unknown sReconDist = %s specified', sNullModel);
        return;
end

% (w - \gammm mP)
mWeighting = mData - gamma * mP;
dist = bmCompareRecon(mMembership1, mMembership2, mImageMat1, mImageMat2);
dist = dist .* mWeighting;


end