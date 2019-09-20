classdef InitRandomImage
    %
    % @author: Jeffrey Chan, 2013
    %

    properties
        m_bHard = true;
    end
    
    methods
        
        function obj = InitRandomImage(bHard)
            obj.m_bHard = bHard;
        end
    
        function [mImage] = initImage(obj, mAdj, mMembership)
        %
        % Initialise the image matrix.
        %
        % Random.
        %
            if obj.m_bHard
                mImage = randi(2, size(mMembership,2)) - 1;
            else
                mImage = rand(size(mMembership,2));
            end
        end % end of function
    end % end of methods
    
    
    
end % end of class


