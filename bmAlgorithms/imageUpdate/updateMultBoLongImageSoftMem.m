
function [mImage] = updateMultBoLongImageSoftMem(mAdj, mImage, mMembership, fDistanceFunc)
%
% Updates mImage according to KKT condition derived multiplicative update
% rule.
%
% Bo Long formulation for Soft cluster membership
%
%
% @author: Jeffrey Chan, 2014
%
    mMemSq = mMembership' * mMembership;

%     mDenom = (mMemSq * mImage * mMemSq);
    
    mImage = mImage .* ((mMembership' * mAdj * mMembership + eps) ./ (mMemSq * mImage * mMemSq + eps));
    
%     mNomin = (mMembership' * mAdj * mMembership);
%     % avoid divide by zero errors
%     mNonZero = mDenom ~= 0;
%     mImage(mNonZero) = mImage(mNonZero) .* (mNomin(mNonZero) ./ mDenom(mNonZero));
%     mImage(~mNonZero) = 0;    
    
%     mImage = mImage .* ((mMembership' * mAdj * mMembership) ./ (mMemSq * mImage * mMemSq));

end % end of function