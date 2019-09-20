
function [mImage] = updateMultBoLongImageSoftMemIDiv(mAdj, mImage, mMembership, fDistanceFunc)
%
% Updates mImage according to KKT condition derived multiplicative update
% rule.
%
% Bo Long formulation for Soft cluster membership
%
%
% @author: Jeffrey Chan, 2014
%
    mImage = mImage .* ((mMembership' * (mAdj ./ (mMembership * mImage * mMembership')) * mMembership) ./ (mMembership' * ones(size(mMembeship,1), size(mMembership,2)) * mMembership));

end % end of function