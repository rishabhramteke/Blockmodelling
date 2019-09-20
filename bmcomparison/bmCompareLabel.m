function [dist] = bmCompareLabel(mMembership1, mMembership2, mImageMat1, mImageMat2)
%
% mMembership1 - matrix of memberships of each vertex
% mMembership2 - matrix of memberships of each vertex
% mImageMat1 - image matrix
% mImageMat2 - image matrix
%
% Compares the blockmodels represented by cPosition1 and cPosition2.
%

% mOverlap = mMembership1' * mMembership2;
% 
% dist = mOverlap * mOverlap * abs(mImageMat1 - mImageMat2);
% % need to normalise by n^2 b/c we didn't do that for mOverlap
% dist = dist / size(mMembership1, 1);

assert(size(mMembership1,1) == size(mMembership2,1));
dist = bmCompareRecon(mMembership1, mMembership2, mImageMat1, mImageMat2) / (size(mMembership1,1)^2);


end % end of function