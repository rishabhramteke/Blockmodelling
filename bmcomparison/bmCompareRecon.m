function [dist] = bmCompareRecon(mMembership1, mMembership2, mImageMat1, mImageMat2, varargin)
%
% Compares the blockmodels represented by cPosition1 and cPosition2,
% reconstruction error.
%
% INPUT:
% mMembership1 - matrix of memberships of each vertex
% mMembership2 - matrix of memberships of each vertex
% mImageMat1 - image matrix
% mImageMat2 - image matrix
% varargin{1} - weighting matrix (adjusted for randomness)
%
% OUTPUT:
% dist          - distance between the blockmodels.
%
% If using this code, please kindly cite the following paper:
% J. Chan, X. Nguyen, W. Liu, C. Leckie, J. Bailey, K. Ramamohanarao and J. Pei. 
% "Structure-aware Distance Measures for Comparing Clusterings in Graphs". 
% In Proceedings of the 18th Pacific-Asia Conference on Knowledge Discovery and Data Mining, May 2014. 
%
% @author: Jeffrey Chan, 2013
%

    assert(size(mMembership1,1) == size(mMembership2,1));

    if length(varargin) == 0
        dist = sum(sum(abs(mMembership1 * mImageMat1 * mMembership1' - mMembership2 * mImageMat2 * mMembership2'))) / (size(mMembership1,1)^2);
    else
        mWeight = varargin{1};
        dist = sum(sum(abs((mMembership1 * mImageMat1 * mMembership1' - mMembership2 * mImageMat2 * mMembership2') .* mWeight))) / (size(mMembership1,1)^2);
    end

end % end of function