function [mMembership] = updateSoftMembershipMult(mAdj, mImage, mMembership)
%
% Updates the soft membership, using multiplicatinve update rule. (Bo Long)
%
%

    mXS = mMembership * mImage;
    mXSt = mMembership * mImage';
    mMembership = mMembership * ( (mAdj' * mXS + mAdj * mXSt) ./ (mXS * (mMembership') * mXSt + mXSt * (mMembership') * mXS) )^0.25; 


end % end of function