function [dist] = bmCompareReconKL(mMembership1, mMembership2, mImageMat1, mImageMat2)
%
% Compares the blockmodels represented by cPosition1 and cPosition2,
% KL reconstruction error.
%
% INPUT:
% mMembership1 - matrix of memberships of each vertex
% mMembership2 - matrix of memberships of each vertex
% mImageMat1 - image matrix
% mImageMat2 - image matrix
%
% OUTPUT:
% dist       - distance between the blockmodels.
%
% If using this code, please kindly cite the following paper:
% J. Chan, X. Nguyen, W. Liu, C. Leckie, J. Bailey, K. Ramamohanarao and J. Pei. 
% "Structure-aware Distance Measures for Comparing Clusterings in Graphs". 
% In Proceedings of the 18th Pacific-Asia Conference on Knowledge Discovery and Data Mining, May 2014. 
%
% @author: Jeffrey Chan, 2013
%

    mAdj1 = mMembership1 * mImageMat1 * mMembership1';
    mAdj2 = mMembership2 * mImageMat2 * mMembership2';

    dist = sum(sum(((mAdj1+eps) .* log2((mAdj1+eps) ./ (mAdj2+eps))) - mAdj1 + mAdj2));

end % end of function


