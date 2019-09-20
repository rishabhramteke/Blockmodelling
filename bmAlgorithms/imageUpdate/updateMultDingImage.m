function [mImage] = updateMultDingImage(mAdj, mImage, mMembership, fDistanceFunc)
%
% Updates the image matrix according to the multiplicative update rule described
% in "Community discovery using nonnegative matrix factorization"
%
%
% @author: Jeffrey Chan, 2013
%

    mXtX = mMembership' * mMembership;
%     mDenom = (mXtX * mImage * mXtX);
%     mNomin = (mMembership' * mAdj * mMembership);
%     % avoid divide by zero errors
%     mNonZero = mDenom ~= 0;
%     mImage(mNonZero) = mImage(mNonZero) .* (mNomin(mNonZero) ./ mDenom(mNonZero));
%     mImage(~mNonZero) = 0;
    
%     mFactor = (mMembership' * mAdj * mMembership) ./ (mXtX * mImage * mXtX);
%     % remove any divide by zero errors
%     mFactor(isnan(mFactor)) = 0;
%     mFactor(isinf(mFactor)) = 0;
%     % update image
    mImage = mImage .* (mMembership' * mAdj * mMembership + eps) ./ (mXtX * mImage * mXtX + eps);
end % end of function