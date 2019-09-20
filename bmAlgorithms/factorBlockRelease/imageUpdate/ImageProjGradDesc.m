classdef ImageProjGradDesc
    %
    % Gradient descent method for updating the image matrix, adjusted with
    % regularisation term.
    %
    % If using this code, please kindly cite the following paper:
    % J. Chan, W. Liu, C. Leckie, J. Bailey and K. Ramamohanarao. 
    % "Discovering Latent Blockmodels in Sparse and Noisy Graphs using Non-Negative Matrix Factorisation."
    % In Proceedings of 22nd ACM International Conference on Information and Knowledge Management, October 2013.            
    %
    % @author: Jeffrey Chan, 2014
    %
    
    properties
        % step size weightings
        m_alpha;
        % calculate alpha each time
        m_bCalAlpha;
        m_beta;
        m_sigma;
    end
   
    methods
    
        function obj = ImageProjGradDesc(varargin)
            
            % parse arguments
            inParser = inputParser;

            % add parameters and set default valuse
            addParameter(inParser, 'initStepSize', 1);
            addParameter(inParser, 'calcAlpha', true);
            addParameter(inParser, 'stepSizeChangeMult', 0.1);
            addParameter(inParser, 'stepStopThres', 0.01);

            
            parse(inParser, varargin{:});
    
            obj.m_alpha = inParser.Results.initStepSize;
            obj.m_bCalAlpha = inParser.Results.calcAlpha;
            obj.m_beta = inParser.Results.stepSizeChangeMult;
            obj.m_sigma = inParser.Results.stepStopThres;    
         
        end
        
        
        function [mImage] = updateImage(obj, mAdj, mImage, mMembership, fObjective)
            %
            % Updates mImage according to a projected gradient descent update rule.
            %
            
            if obj.m_bCalAlpha
                alpha = (size(mImage,1) / size(mAdj,1))^2;
            else
                alpha = obj.m_alpha;
            end
            
            beta = obj.m_beta;
            sigma = obj.m_sigma;

            mGrad = fObjective.gradient(mAdj, mImage, mMembership);
            mImageNew = obj.posProjection(mImage + alpha * mGrad);
            currDis = fObjective.distance(mAdj, mImage, mMembership);

            if fObjective.distance(mAdj, mImageNew, mMembership) - currDis > sigma * sum(sum(mGrad .* (mImageNew - mImage)))
                % compute new alpha and new image and gradient
                alpha = alpha * beta;
                mImageNew = obj.posProjection(mImage + alpha * mGrad);
                % while difference in magnitude isn't a "minimum", decrease alpha
                while fObjective.distance(mAdj, mImageNew, mMembership) - currDis > sigma * sum(sum(mGrad .* (mImageNew - mImage)))
                    % compute new alpha and new image and gradient
                    alpha = alpha * beta;
                    mImageNew = obj.posProjection(mImage + alpha * mGrad);
                end
            else
                % compute new alpha and new image and gradient
                alpha = alpha / beta;
                imageDiffLimit = 10^-9;
                mImageNew = obj.posProjection(mImage + alpha * mGrad);
                while abs(sum(sum((mImageNew - mImage)))) > imageDiffLimit && (fObjective.distance(mAdj, mImageNew, mMembership) - currDis <= sigma * sum(sum(mGrad .* (mImageNew - mImage))))
                    if (isinf(alpha))
                        error('alpha is infinity');
                    end
                     % compute new alpha and new image and gradient
                    alpha = alpha / beta;
                    mImageNew = obj.posProjection(mImage + alpha * mGrad);
                end
            end % end of if


        end % end of function
    
        

        function mMatNew = posProjection(~, mMat)
        %
        % Projects each entry of mMat to non-negative values.
        %
            mMatNew = mMat;
            mMatNew(mMat < 0) = 0;
        end        
        
    end % end of methods
    
    
end % end of class



