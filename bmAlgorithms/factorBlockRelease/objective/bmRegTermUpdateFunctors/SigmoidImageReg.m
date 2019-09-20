classdef SigmoidImageReg < handle
    %
    % Sigmoid image regularisation function object.
    %
    % If using this code, please kindly cite the following paper:
    % J. Chan, W. Liu, C. Leckie, J. Bailey and K. Ramamohanarao. 
    % "Discovering Latent Blockmodels in Sparse and Noisy Graphs using Non-Negative Matrix Factorisation."
    % In Proceedings of 22nd ACM International Conference on Information and Knowledge Management, October 2013.   
    %
    %
    % @author: Jeffrey Chan, 2013
    %    
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        m_gamma = 0.25;
        m_beta_mu = 500;
        m_mTau = 0.5;
    end % end of properties
    
    
    methods
            
        function obj = SigmoidImageReg(varargin)
            if ~isempty(varargin)
                if length(varargin) >= 1
                    obj.m_gamma = varargin{1};
                end
                if length(varargin) >= 2
                    obj.m_beta_mu = varargin{2};
                end
                if length(varargin) == 3
                    obj.m_mTau = varargin{3};
                end 
            end
        end % end of function
        
        
        function adjustIdealShift(obj, mTau)
            %
            % Adjust the tau component of the sigmoid function.
            %
            
            obj.m_mTau = mTau;
        end % end of function
        
        
        function dis = computeImageRegObj(obj, mImage)
            %
            % Compute the contribution of the regularisation term to the overall
            % objective.
            %
            
            mIdeal = 1 ./ (1 + obj.m_gamma * exp(-obj.m_beta_mu * (mImage - obj.m_mTau)));
            mApprox = mIdeal - mImage;
            dis = trace(mApprox * mApprox');   
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