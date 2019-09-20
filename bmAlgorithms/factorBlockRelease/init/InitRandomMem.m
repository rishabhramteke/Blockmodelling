classdef InitRandomMem
    %
    % Randomly initialise membership matrix.
    % 
    % If using this code, please kindly cite the following paper:
    % J. Chan, W. Liu, C. Leckie, J. Bailey and K. Ramamohanarao. 
    % "Discovering Latent Blockmodels in Sparse and Noisy Graphs using Non-Negative Matrix Factorisation."
    % In Proceedings of 22nd ACM International Conference on Information and Knowledge Management, October 2013.
    %
    % @author: Jeffrey Chan, 2014
    %

    properties
        % hard cluster initialisation or not
        m_bHard = true;
    end
    
    methods
        
        function obj = InitRandomMem(bHard)
            obj.m_bHard = bHard;
        end
    
        function [mMembership] = initMembership(obj, mAdj, k)
        %
        % Initialise the membership matrix.
        %
        % Random.
        %
            n = size(mAdj,1);
        
            if obj.m_bHard
                mMembership = randomInitMembership(n,k);
            else
                mMembership = rand(n,k);
            end
        end % end of function
    end % end of methods
    
    
end % end of class


function [mMembership] = randomInitMembership(n, k)
    % initialise mMembership
    mMembership = zeros(n, k);


    % vRoleNum = (randfixedsum(k, 1, n, 1, n))';
    
    assert(n >= k);
    
    % assign at least one vertex to each position
    vSample = randsample(n, k);
    vRemaining = logical(true(1,n));
    vRemaining(vSample) = false;
    
    
    for v = 1 : size(vSample,1)
        mMembership(vSample(v), v) = 1;
    end
    

    vRemainIndices = find(vRemaining);
    vNewPos = randi(k, length(vRemainIndices), 1);
    mIndices = logical(sparse(vRemainIndices, vNewPos, ones(length(vRemainIndices),1), size(mMembership,1), size(mMembership,2)));
    mMembership(mIndices) = 1;


    % convert to sparse representation
    mMembership = sparse(mMembership);
    
    for v = 1 : n
        if isempty(find(mMembership(v,:) > 0))
            fprintf('vertex %d missing assignment', v);
        end
    end

end        

