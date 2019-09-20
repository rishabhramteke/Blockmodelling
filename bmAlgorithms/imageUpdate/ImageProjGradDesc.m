classdef ImageProjGradDesc
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
        m_maxValue;
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
            addParameter(inParser, 'maxValue', 0.01);

            
            parse(inParser, varargin{:});
    
            obj.m_alpha = inParser.Results.initStepSize;
            obj.m_bCalAlpha = inParser.Results.calcAlpha;
            obj.m_beta = inParser.Results.stepSizeChangeMult;
            obj.m_sigma = inParser.Results.stepStopThres;    
            obj.m_maxValue = inParser.Results.maxValue;    
         
        end
        
        
        function [mImageNew] = updateImage(obj, mAdj, mImage, mMembership, fObjective)
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

            tinyAlpha = 1^-5;
            
        %     if fObjective.distance(mAdj, mImageNew, mMembership) - currDis > sigma * mGrad' * (mImageNew - mImage)
            if fObjective.distance(mAdj, mImageNew, mMembership) - currDis > sigma * sum(sum(mGrad .* (mImageNew - mImage)))
                % compute new alpha and new image and gradient
                alpha = alpha * beta;
                mImageNew = obj.posProjection(mImage + alpha * mGrad);
                % while difference in magnitude isn't a "minimum", decrease alpha
                while alpha > tinyAlpha && fObjective.distance(mAdj, mImageNew, mMembership) - currDis > sigma * sum(sum(mGrad .* (mImageNew - mImage))) 
        %             disp(sprintf('alpha = %f', alpha));
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
        %             disp(sprintf('alpha = %f', alpha));
                     % compute new alpha and new image and gradient
                    alpha2 = alpha / beta;
                    if (isinf(alpha2))
%                         error('alpha is infinity');
                        mImageNew = obj.posProjection(mImage + alpha * mGrad);
                        break;
                    end    
                    
                    mImageNew = obj.posProjection(mImage + alpha2 * mGrad);
                    
                    alpha = alpha2;
                end
                
            end % end of if


        end % end of function
    
        


        function mMatNew = posProjection(obj, mMat)
        %
        % Projects each entry of mMat to non-negative values.
        %
            mMatNew = mMat;
            mMatNew(mMat < 0) = 0;
            mMatNew(isnan(mMat)) = obj.m_maxValue;
            mMatNew(mMat > obj.m_maxValue) = obj.m_maxValue;
        end        
        
    end % end of methods
    
    
end % end of class



