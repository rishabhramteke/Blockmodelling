classdef ImagePurtity
    %
    % Comput the purity of image matrix (i.e., how far from 0/1 it is).
    %
    % @author: Jeffrey Chan, 2014
    %
    
    properties
    end
    
    
    methods
        function [obj] = ImagePurtity()

        end
        
        
        function [dis] = distance(~, mImage, ~)
            mVal = zeros(size(mImage,1), size(mImage,2));
            mNonZero = mImage ~= 0;
            % compute entropy for nonzero entries
            mVal(mNonZero) = -mImage(mNonZero) .* log2(mImage(mNonZero));
            
            dis = sum(sum(mVal));
        end
        

    end
    
    
end % end of class