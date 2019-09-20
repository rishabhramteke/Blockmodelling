function [denDis] = bmUnweightedRecon(mMembership1, mImage1, mData1, mMembership2, mImage2, mData2,...
    sReconDist)
%
% Computes the density distances between two blockmodels.
%
% INPUT:
% mMembership1  - matrix of memberships of each vertex
% mImageMat1    - image matrix 
% mData1        - adjacency matrix
% mMembership2  - matrix of memberships of each vertex
% mImageMat2    - image matrix
% mData1        - adjacency matrix
% sReconDist    - Reconstruction distance to use.
%
% OUTPUT:
% denDist       - distance between the blockmodels.
%
% If using this code, please kindly cite the following paper:
% J. Chan, X. Nguyen, W. Liu, C. Leckie, J. Bailey, K. Ramamohanarao and J. Pei. 
% "Structure-aware Distance Measures for Comparing Clusterings in Graphs". 
% In Proceedings of the 18th Pacific-Asia Conference on Knowledge Discovery and Data Mining, May 2014. 
%
% @author: Jeffrey Chan, 2013
%



switch sReconDist
    case 'edgeRecon'
        % compute density distance
        denDis = bmCompareRecon(mMembership1, mMembership2, mImage1, mImage2);
    case 'blockRecon'
        % block distance (same as edge densit distance
        denDis = bmCompareRecon(mMembership1, mMembership2, mImage1, mImage2);
    case 'edit'
        % compute edit distance
        mXor = xor(mData1, mData2);
        denDis = size(find(mXor),1);
    case 'edgeReconKL'
        denDis = bmCompareReconKL(mMembership1, mMembership2, mImage1, mImage2);
    otherwise
        warning('Unknown sReconDist = %s specified', sReconDist);
        return;
end


end % end of function




