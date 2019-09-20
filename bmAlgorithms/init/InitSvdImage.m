classdef InitSvdImage
    %
    % Initialisation of the image matrix based on the SVD method described in
    % 'SVD based initialization: A head start for nonnegative matrix
    % factorization." 
    %
    % @author: Jeffrey Chan, 2014
    %

    properties
        m_sDenseOption = 'None';
    end
    
    methods
        
        function obj = InitSvdImage(varargin)
            % parse arguments
            inParser = inputParser;

            % add parameters and set default valuse
            addParameter(inParser, 'denseOption', 'Da');

            parse(inParser, varargin{:});
    
            obj.m_sDenseOption = inParser.Results.denseOption;           
        end % end of function
    
        
        function [mImage] = initImage(obj, mAdj, mMembership)
        %
        % Initialise the image matrix using nonnegative projects of SVD
        % decomposition.
        %
            k = size(mMembership,2);
            [mU, mS, mV] = svds(mAdj, k);
            
        
            if obj.m_bHard
                mImage = randi(2, size(mMembership,2)) - 1;
            else
                mImage = rand(size(mMembership,2));
            end
        end % end of function
        
    end % end of methods
    
    
    
end % end of class


