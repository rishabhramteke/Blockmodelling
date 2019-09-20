classdef ImageMultHardSumConstraintRepara
    %
    % Multplicative soft image update function with sum(sum(mImage)) = m*k^2/n^2 
    % Uses the reparameterisation approach.
    %
    % Introduce well as L1 and L2 regularisation as options.
    %
    % @author: Jeffrey Chan, 2014
    %    
    
    properties
       % regularisation parameters
       m_paraL1MembershipReg = 0;
       m_paraL2MembershipReg = 0;
    end
    
    
    methods
        
        
        function obj = ImageMultHardSumConstraintRepara(varargin)
        %
        % INPUT:
        % varargin{1}:
        % varargin{2}:
        %
              
            if length(varargin) >= 1
                assert(varargin{1} >= 0);
                obj.m_paraL1MembershipReg = varargin{1};
            end
            if length(varargin) == 2
                assert(varargin{2} >= 0);
                obj.m_paraL2MembershipReg = varargin{2};
            end
        end % end of function
        
        
        function [mImage] = updateImage(obj, mAdj, mImage, mMembership, fDistanceFunc)
        %
        % Using multiplicatinve update rule.
        %
        %

            mXtX = mMembership' * mMembership;
            
            % regularisation gradient
            mRegGraD = fDistanceFunc.imageRegGradient(mImage);
            % update image
            mPos = mMembership' * mAdj * mMembership + fDistanceFunc.m_regPara * mRegGraD + obj.m_paraL1MembershipReg + obj.m_paraL2MembershipReg * mImage;
            mNeg = mXtX * mImage * mXtX;
            imageSum = sum(sum(mImage));
            sum1 = sum(sum(mPos .* mImage));
            sum2 = sum(sum(mNeg .* mImage));
            
            mImage = mImage .* ((mNeg * imageSum + sum1) ./ (mPos * imageSum + sum2 + eps));
            
            
            % normalise
            vertNum = size(mMembership,1);
            posNum = size(mMembership,2);
            edgeNum = nnz(mMembership);
            if edgeNum ~= 0
                mImage = mImage ./ (sum(sum(mImage)) * (vertNum^2 / (posNum^2 * edgeNum)));
            else
                mImage = zeros(sizeof(mImage,1), sizeof(mImage,2));
            end
            
            mImage = max(eps, mImage);
            
        end % end of function
        
        
    end % end of methods
    
end