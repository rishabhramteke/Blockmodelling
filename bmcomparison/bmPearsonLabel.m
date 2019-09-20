function dist = bmPearsonLabel(mMembership1, mMembership2, mImage1, mImage2, mAdj1, mAdj2)
%
% Computes the pearson similarity proposed in waserman's book.
%


% compute the mean of the two matrices
mean1 = mean(mean(mAdj1));
mean2 = mean(mean(mAdj2));

% assume partitional
mApproxAdj1 = mMembership1 * mImage1 * mMembership1';
mApproxAdj2 = mMembership2 * mImage2 * mMembership2';

% compute the cross and variances
mFirstMom1 = mApproxAdj1 - mean1;
mFirstMom2 = mApproxAdj2 - mean2;

crossVar = sum(sum(mFirstMom1 .* mFirstMom2));
var1 = norm(mFirstMom1, 'fro');
var2 = norm(mFirstMom2, 'fro');


% check if any of the first moment matrices are 0
% if so, dist = 0
if (var1 == 0 || var2 == 0)
    dist = 0;
else
    dist = crossVar / (var1 * var2);
end

end