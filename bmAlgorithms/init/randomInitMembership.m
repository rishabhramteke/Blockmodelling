function [mMembership] = randomInitMembership(n, k)
%
% Randomly generate a valid (hard/partitional) membership matrix
%

    % initialise mMembership
    mMembership = zeros(n, k);


    % vRoleNum = (randfixedsum(k, 1, n, 1, n))';
    
    assert(n >= k);
    
    % assign at least one vertex to each position
    vSample = randsample(n, k);
    vRemaining = logical(true(1,n));
    vRemaining(vSample) = false;
    
%     vRemaining = setdiff([1:n], vSample);
    
    for v = 1 : size(vSample,1)
        mMembership(vSample(v), v) = 1;
    end
    

    vRemainIndices = find(vRemaining);
    vNewPos = randi(k, length(vRemainIndices), 1);
    mIndices = logical(sparse(vRemainIndices, vNewPos, ones(length(vRemainIndices),1), size(mMembership,1), size(mMembership,2)));
    mMembership(mIndices) = 1;
    
%     for v = 1 : length(vRemainIndices)
%         newK = randsample(k, 1); 
%         mMembership(vRemainIndices(v), newK) = 1;
%     end



    % convert to sparse representation
    mMembership = sparse(mMembership);
    
    for v = 1 : n
        if isempty(find(mMembership(v,:) > 0))
            disp(sprintf('vertex %d missing assignment', v));
        end
    end

end