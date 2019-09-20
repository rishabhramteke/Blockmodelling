classdef GaussianImageReg
    %
    % Gaussian image regularisation function object.
    %
    %
    % @author: Jeffrey Chan, 2013
    %    
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        % default values (make sure they match with EucTriFactBM)
        m_
    end % end of properties
    
    
    methods
            
        function obj = GaussianImageReg(varargin)

        end % end of function
        
        
        function dis = computeImageRegObj(obj, mImage)
            %
            % Compute the contribution of the regularisation term to the overall
            % objective.
            %
            
            % (exp(-0.5 * delta))^2 = exp(-delta)
            mIdeal = exp(-(mImage - obj.m_mGaussMean).^2);
            dis = sum(sum(mIdeal));
%             dis = trace(mIdeal * mIdeal');   
        end % end of function
        
        
        function [mImageRegGrad] = computeImageRegGrad(obj, mImage)
            %
            % Compute the gradient of the image regularisation term.
            % Does not factor in the regularisation weight.
            %
            
            mExp = obj.m_gamma * exp(-obj.m_beta_mu *(mImage - obj.m_mTau));
            mExpSq = (1 + mExp) .* (1 + mExp);
            mThird = 1 + mExp .* (1 + obj.m_beta_mu * mImage) ./ mExpSq;
            
            mExpCb = mExpSq .* (1 + mExp);
            mFifth = 2 * obj.m_beta_mu * (mExp ./ mExpCb);
            
            mImageRegGrad = mThird + 2 * mImage + mFifth;
        end % end of function()
    
        
        
        function dis = computeImageRegCoord(obj, mImage, stepSize, basisRow, basisCol)
            %
            % Compute the image regularisation objective for coordinate descent.
            %
            
            disSq = stepSize * stepSize;
                
            mX2 = mImage - obj.m_mTau;
            % compute the term where we don't add t
            mExp = 1 ./ (1 + obj.m_gamma * exp(-obj.m_beta_mu * mX2));
            % subtract the exponential term that isn't included in the sum
            mExp(basisCol, basisRow) = 0;
            
            % compute the exp term where we add stepSize
            singleExp = 1 ./ (1 + obj.m_gamma * exp(-obj.m_beta_mu * (stepSize + mX2(basisRow, basisCol))));
            
            dis = (disSq + singleExp * singleExp - 2 * stepSize * singleExp + 2 * stepSize * mImage(basisRow, basisCol)...
                - 2 * mImage(basisRow, basisCol) * singleExp + trace(mExp' * (mExp + mImage)) + trace(mImage' * mImage));
        end % end of function
        
    end % end of methods
    
end % end of class