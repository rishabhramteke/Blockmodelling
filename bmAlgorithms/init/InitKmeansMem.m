classdef InitKmeansMem
    %
    % Spherical K-means based initialisation.
    %
    % @author: Jeffrey Chan, 2014
    %

    properties
        m_bHard = true;
    end
    
    methods
        
        function obj = InitKmeansMem(bHard)
            obj.m_bHard = bHard;
        end
    
        function [mMembership] = initMembership(obj, mAdj, k)
        %
        % Initialise the membership matrix.
        %
            n = size(mAdj,1);
            mMembership = zeros(n, k);

        end % end of function
    end % end of methods
    
    
    
end % end of class




