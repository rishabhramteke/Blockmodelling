classdef EucTriFactBMAdj < EucTriFactBM
    % 
    % Euclidean tri-factorisation with sigmiod constraint and has adjusted
    % edge mismatch penalties.
    %
    %
    % @author: Jeffrey Chan, 2013
    %
    
    
    properties 
        % pre-computed weight matrix.
        m_mWeight;
        m_vWeight;
    end % end of properties

    
    
    methods
        function obj = EucTriFactBMAdj(mWeight, varargin)
            %
            % constructor
            % mWeight - weight adjustment matrix.
            %
            
            % call super class
            obj@EucTriFactBM(varargin{:});
           
            obj.m_mWeight = mWeight;
            obj.m_vWeight = mWeight(:);
        end
        
        
        function [dis] = distance(obj, mAdj, mImage, mMembership)
            %
            % Compute the approximation distance between mAdj and the
            % approximated blockmodel + regularisation term.
            %

            dis = distance@EucTriFactBM(obj, mAdj, mImage, mMembership, obj.m_mWeight);
        end % end of distance()
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function [mGrad] = gradient(obj, mAdj, mImage, mMembership)
            %
            % Compute the gradient of approximation function.
            %
            
            % rescaling factor for regularisation term
%             rescalFactor = (size(mAdj,1)^2 / size(mImage,1)^2);
            
            % this is the term that is different from EucTriFactBM
            mTerm = 2 * mMembership' * ((mMembership * mImage * mMembership' - mAdj) .* obj.m_mWeight.^2) * mMembership;

            
            mRegGrad = obj.imageRegGradient(mImage);
               
            % add all the pieces intogether
%             mGrad = mTerm + obj.m_regPara * rescalFactor * mRegGrad;
            mGrad = mTerm + obj.m_regPara * mRegGrad;
        end % end of function gradient()
        
        
        function vX0 = computeApprox(obj, mAdj, mImage, mMembership)
            vX0 = obj.m_fComputeCoordDesc.computeWeightedApprox(mAdj, mImage, mMembership, obj.m_mWeight);
        end % end of function
        
        

        function [ vX1 ] = computeX1(obj, mMembership, basisRow, basisCol)
           vX1 =  obj.m_fComputeCoordDesc.computeWeightedX1(mMembership, basisRow, basisCol, obj.m_mWeight);
        end
        
        
        function vX0 = updateX0Add(obj, vX0, mX1, stepSize)
            vX0 = obj.m_fComputeCoordDesc.updateWeightedX0Add(vX0, mX1, stepSize, obj.m_vWeight);
        end                
        
        
        
        
        function [dis] = coordDescDistance(obj, mAdj, mImage, mMembership, stepSize, basisRow, basisCol)
            %
            % Compute the objective value of the coordinate descent expression.
            % TODO: we can reuse most of the calculations again, need to figure
            % out how to store them between calls.
            %
            
            % rescaling factor for regularisation term
%             rescalFactor = (size(mAdj,1)^2 / size(mImage,1)^2);
            
            % construct the unit basis
%             mUnitBasis = zeros(size(mImage,1), size(mImage,2));
%             mUnitBasis(basisRow, basisCol) = 1;
%             mUnitBasis = sparse(mUnitBasis);
%             
%             mX0 = (mAdj - mMembership * mImage * mMembership') .* obj.m_mWeight;
%             mX1 = (mMembership * mUnitBasis * mMembership') .* obj.m_mWeight;
%             disSq = stepSize * stepSize;
%             % Tr(X0' X0)
%             dis = trace(mX0' * mX0);
%             % Tr(X0' X1)
%             dis = dis - 2 * stepSize * trace(mX0' * mX1);
%             % Tr(X1' X1)
%             dis = dis + disSq * trace(mX1' * mX1);
            
           
            dis = coordDescObj(mAdj, mImage, mMembership, stepSize, basisRow, basisCol, obj.m_mWeight);            
            
%             mX2 = mImage - obj.m_mTau;
%             % compute the term where we don't add t
%             mExp = 1 ./ (1 + obj.m_gamma * exp(-obj.m_beta_mu * mX2));
%             % subtract the exponential term that isn't included in the sum
%             mExp(basisCol, basisRow) = 0;
%             
%             % compute the exp term where we add stepSize
%             singleExp = 1 ./ (1 + obj.m_gamma * exp(-obj.m_beta_mu * (stepSize + mX2(basisRow, basisCol))));
            
%             % compute the whole sum
%             dis = dis + obj.m_regPara * rescalFactor * (disSq + singleExp * singleExp - 2 * stepSize * singleExp + 2 * stepSize * mImage(basisRow, basisCol)...
%                 - 2 * mImage(basisRow, basisCol) * singleExp + trace(mExp' * (mExp + mImage)) + trace(mImage' * mImage));            
%             dis = dis + obj.m_regPara * rescalFactor * obj.imageRegCoordDesc(mImage, stepSize, basisRow, basisCol);            
            dis = dis + obj.m_regPara * obj.imageRegCoordDesc(mImage, stepSize, basisRow, basisCol);            
        end % end of function coordDescDistance
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        
        
        function [mGrad] = posGradient(obj, mAdj, mImage, mMembership)
            %
            % Compute the gradient of the position matrix.
            %
            
            mApproxDiff = (mMembership * mImage * mMembership' - mAdj) .* obj.m_mWeight .* obj.m_mWeight;
            
            mGrad = 2* (mApproxDiff' * mMembership * mImage + mApproxDiff * mMembership * mImage');
            
            mGradExtra = obj.posGradientExtraTerm(mAdj, mImage, mMembership);
            if ~isnan(mGradExtra)
                mGrad = mGrad + mGradExtra;
            end
        end % end of function posGradient()        
        
    
    end




end % end of class