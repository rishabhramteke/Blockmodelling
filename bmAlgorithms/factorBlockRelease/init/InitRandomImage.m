classdef InitRandomImage
    %
    % Randomly initialise image matrix.
    %
    % If using this code, please kindly cite the following paper:
    % J. Chan, W. Liu, C. Leckie, J. Bailey and K. Ramamohanarao. 
    % "Discovering Latent Blockmodels in Sparse and Noisy Graphs using Non-Negative Matrix Factorisation."
    % In Proceedings of 22nd ACM International Conference on Information and Knowledge Management, October 2013.    
    %
    % @author: Jeffrey Chan, 2013
    %

    properties
        % hard cluster initialisation or not
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


