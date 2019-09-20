function [posDis, denDis] = bmOverlapDist(mMembership1, mData1, mMembership2, mData2,...
    sReconDist, sPosDist)
%
% Computes the positional and density distances between clustering/bm 1 and
% 2.
%
% sReconDist - Reconstruction distance to use.
% sPosDist - Position distance to use.
%


denDis = bmOverlapDenDist(mMembership1, mData1, mMembership2, mData2, sReconDist);

posDis = bmOverlapPosDist(mMembership1, mMembership2, sPosDist);



end % end of function




