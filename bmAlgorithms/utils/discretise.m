function [mNew] = discretise(mMembership)
    %
    % Discretise the input matrix.
    % Assign 1 to element of each row that is max.
    %
    % @author: Jeffrey Chan, 2013
    %
    
    [~, vMaxIndex] = max(mMembership, [], 2);
    mB = cat(2, [1:size(mMembership,1)]', vMaxIndex);
    mIndices = logical(sparse(mB(:,1), mB(:,2), 1, size(mMembership,1), size(mMembership,2)));
    mNew = double(mIndices);
end
