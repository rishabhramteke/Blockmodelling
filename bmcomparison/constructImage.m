function [mImage] = constructImage(mMembership, mData)
%
% Constructs an image matrix for mMembership and mData.
%

mImage = zeros(size(mMembership,2), size(mMembership,2));
mBlockSize = zeros(size(mMembership,2));
% get all nonzero entries in MData1
[R, C, Val] = find(mData);
assert(size(R,1) == size(C,1));
for v = 1 : size(R,1)
    % look up position for row and col vertices
    vrowPos = find(mMembership(R(v),:) > 0);
    vcolPos = find(mMembership(C(v),:) > 0);
    % TODO: assume partitional atm
    
    % TODO: assume hard partitioning, by adding 1 to densities
    mImage(vrowPos(1), vcolPos(1)) = mImage(vrowPos(1), vcolPos(1)) + 1;
end

% compute the size of each block
for k = 1 : size(mMembership,2)
    mBlockSize(k) = size(find(mMembership(:,k) > 0),1);
end

% normalise image matrix with block sizes
for r = 1 : size(mImage,1)
    for c = 1 : size(mImage,2)
        mImage(r,c) = mImage(r,c) / (mBlockSize(r) * mBlockSize(c));
    end
end


end