classdef ImageMultHardSumConstraint
    %
    % Multplicative soft image update function with sum(sum(mImage)) = m*k^2/n^2 
    % Uses the Lagragian approximation approach.
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
        
        
        function obj = ImageMultHardSumConstraint(varargin)
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
            
            mNegPosRatio = mNeg ./ (mPos + eps);
            mF = sum(sum(mImage ./ (mPos + eps)));
            mJ = sum(sum(mImage .* mNegPosRatio));
            
            vertNum = size(mMembership,1);
            posNum = size(mMembership,2);
            edgeNum = nnz(mMembership);            
            targetImageSum = (posNum^2 * edgeNum) / vertNum^2;
            
            % update
            mImage = mImage .* ((mNeg * mF + targetImageSum) ./ (mPos * mF + mJ + eps));
            
            mImage = max(eps, mImage);
            
        end % end of function
        
        
    end % end of methods
    
end