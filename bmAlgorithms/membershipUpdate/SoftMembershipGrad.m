classdef SoftMembershipGrad
    %
    % @author: Jeffrey Chan, 2013
    %
    
    properties
        % object function functor
        m_fObjective;
        % step size weightings
        m_alpha;
        % calculate alpha each time
        m_bCalAlpha;
        m_beta1;
        m_beta2;
        m_sigma1;
        m_sigma2;
    end
   
    methods
    
        function obj = SoftMembershipGrad(fObjective, varargin)
            obj.m_fObjective = fObjective;
            
            % parse arguments
            inParser = inputParser;

            % add parameters and set default valuse
            addParameter(inParser, 'initStepSize', 1);
            addParameter(inParser, 'calcAlpha', true);
            addParameter(inParser, 'stepSizeChangeMult1', 0.01);
            addParameter(inParser, 'stepSizeChangeMult2', 0.1);    
            addParameter(inParser, 'stepStopThres1', 0.1);
            addParameter(inParser, 'stepStopThres2', 0.01);
            
            parse(inParser, varargin{:});
    
            obj.m_alpha = inParser.Results.initStepSize;
            obj.m_bCalAlpha = inParser.Results.calcAlpha;
            obj.m_beta1 = inParser.Results.stepSizeChangeMult1;
            obj.m_beta2 = inParser.Results.stepSizeChangeMult2;
            obj.m_sigma1 = inParser.Results.stepStopThres1;    
            obj.m_sigma2 = inParser.Results.stepStopThres2;               
        end
        
        
        function [mMembershipNew] = updateMembership(obj, mAdj, mImage, mMembership)
            %
            % Updates the membership using gradient descent.
            %
            

            % start off roughly in the scale of [0,1] for M
            alpha =  obj.m_alpha;
            beta1 = obj.m_beta1;
            beta2 = obj.m_beta2;
            % larger sigma, then go slower
            sigma1 = obj.m_sigma1;
            sigma2 = obj.m_sigma2;
    
            mGrad =  obj.m_fObjective.posGradient(mAdj, mImage, mMembership);
            
            if obj.m_bCalAlpha
                alpha = 1 / sum(sum(mGrad));
            end            
            
            mMembershipNew = obj.posProjection(mMembership - alpha * mGrad);
            currDis = obj.m_fObjective.distance(mAdj, mImage, mMembership);
            

            
            
            [alpha, mMembershipNew] = obj.findStepSize(alpha, beta1, sigma1, mAdj, mImage, mGrad, mMembership, mMembershipNew, currDis);
            [~, mMembershipNew] = obj.findStepSize(alpha, beta2, sigma2, mAdj, mImage, mGrad, mMembership, mMembershipNew, currDis);
                        
            % this version ensures we don't divide by 0, but slower
%             mMembershipNew = spdiags(spfun(@(x) 1 ./ x, sum(mMembershipNew, 2)), 0, size(mMembershipNew,1), size(mMembershipNew,1)) * mMembershipNew;
        end % end of function
    
        


        function mMatNew = posProjection(~, mMat)
        %
        % Projects each entry of mMat to non-negative values.
        %
            mMatNew = mMat;
            mMatNew(mMat < 0) = 0;
        end        
        
        
        function [alpha, mMembershipNew] = findStepSize(obj, startAlpha, beta, sigma, mAdj, mImage, mGrad, mMembership, mMembershipNew, origDis)
                
            alpha = startAlpha;
    
%             if obj.m_fObjective.distance(mAdj, mImage, mMembershipNew) - currDis > sigma * mGrad' * (mMembershipNew - mMembership)
            if obj.m_fObjective.distance(mAdj, mImage, mMembershipNew) - origDis > sigma * sum(sum(mGrad .* (mMembershipNew - mMembership)))
%                 disp('decrease alpha');
                % compute new alpha and new image and gradient
                alpha = alpha * beta;
                mMembershipNew = obj.posProjection(mMembership - alpha * mGrad);
                % row normalise the rows
%                 mMembershipNew = spdiags(1 ./ sum(mMembershipNew, 2), 0, size(mMembershipNew,1), size(mMembershipNew,1)) * mMembershipNew;                
                % while difference in magnitude isn't a "minimum", decrease alpha
                while obj.m_fObjective.distance(mAdj, mImage, mMembershipNew) - origDis > sigma * sum(sum(mGrad .* (mMembershipNew - mMembership)))
        %             disp(sprintf('alpha = %f', alpha));
                    % compute new alpha and new image and gradient
%                     disp('decrease alpha');
                    alpha = alpha * beta;
                    mMembershipNew = obj.posProjection(mMembership - alpha * mGrad);
                    % row normalise the rows
%                     mMembershipNew = spdiags(1 ./ sum(mMembershipNew, 2), 0, size(mMembershipNew,1), size(mMembershipNew,1)) * mMembershipNew;                    
                end
                % unroll alpha
                alpha = alpha / beta;
            else
%                 disp('increase alpha');
                % compute new alpha and new image and gradient
                alpha = alpha / beta;
%                 membershipDiffLimit = 10^-3;
                mMembershipNew = obj.posProjection(mMembership - alpha * mGrad);
                % row normalise the rows
%                 mMembershipNew = spdiags(1 ./ sum(mMembershipNew, 2), 0, size(mMembershipNew,1), size(mMembershipNew,1)) * mMembershipNew;
%                 while abs(sum(sum((mMembershipNew - mMembership)))) > membershipDiffLimit && (obj.m_fObjective.distance(mAdj, mImage, mMembershipNew) - origDis <= sigma * sum(sum(mGrad .* (mMembershipNew - mMembership))))
                while (obj.m_fObjective.distance(mAdj, mImage, mMembershipNew) - origDis <= sigma * sum(sum(mGrad .* (mMembershipNew - mMembership))))
        %             disp(sprintf('alpha = %f', alpha));
                    if (isinf(alpha))
                        error('alpha has become infinity');
                    end
                    
%                     disp('increase alpha');
                    % compute new alpha and new image and gradient
                    alpha = alpha / beta;
                    mMembershipNew = obj.posProjection(mMembership + alpha * mGrad);
                    % row normalise the rows
%                     mMembershipNew = spdiags(1 ./ sum(mMembershipNew, 2), 0, size(mMembershipNew,1), size(mMembershipNew,1)) * mMembershipNew;
                end
                % unroll last alpha
                alpha = alpha * beta;
            end % end of if
            
        end
        
    end % end of methods
    
    
end % end of class



