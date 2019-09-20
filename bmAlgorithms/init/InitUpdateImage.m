classdef InitUpdateImage
    %
    % @author: Jeffrey Chan, 2013
    %

    properties
        m_fUpdateMethod = NaN;
    end % end of properties
    
    
    methods
        function obj = InitUpdateImage(fUpdateMethod)
            obj.m_fUpdateMethod = fUpdateMethod;
        end
        
        function [mImage] = initImage(obj, mAdj, mMembership)
        %
        % Initialise the image matrix by running the image update function once.
        %

            mImage = obj.m_fUpdateMethod(mAdj, [], mMembership, NaN);
        end % end of function
    end % end of methods

end % end of class definition
