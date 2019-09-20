function [dist] = bmCompareAlign(mMembership1, mMembership2, mImageMat1, mImageMat2)
%
% mMembership1 - matrix of memberships of each vertex
% mMembership2 - matrix of memberships of each vertex
% mImageMat1 - image matrix
% mImageMat2 - image matrix
%
% Compares the blockmodels represented by cPosition1 and cPosition2,
% edge alignment distance.
%

dist= sum(sum(abs(mMembership1 * mImageMat1 * mMembership1' - mMembership2 * mImageMat2 * mMembership2')));



end % end of function