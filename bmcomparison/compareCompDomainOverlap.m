function [posDis, denDis] = compareCompDomainOverlap(mMembership1, mData1, mMembership2, mData2)
%
% Computes the positional and density distances between clustering/bm 1 and
% 2.
%

% construct image matrix (assume partitional for now)
mImage1 = constructImage(mMembership1, mData1);
mImage2 = constructImage(mMembership2, mData2);


% compute density distance
denDis = bmCompareLabel(mMembership1, mMembership2, mImage1, mImage2);
% convert to vector of membership
% TODO: assuming hard partitionoing
vMembership1 = zeros(size(mMembership1,1), 1);
vMembership2 = zeros(size(mMembership2,1), 1);
for v = 1 : size(mMembership1,1)
   vPos = find(mMembership1(v,:) > 0);
   vMembership1(v) = vPos(1);
end
for v = 1 : size(mMembership2,1)
   vPos = find(mMembership2(v,:) > 0);
   vMembership2(v) = vPos(1);
end
posDis = varOfInfo(vMembership1, vMembership2);



end % end of function




