classdef EucTriFactAdj < EucTriFact
    % 
    % Euclidean tri-factorisation frobenius objective that has adjusted
    % edge mismatch penalties.
    %
    % If using this code, please kindly cite the following paper:
    % J. Chan, W. Liu, C. Leckie, J. Bailey and K. Ramamohanarao. 
    % "Discovering Latent Blockmodels in Sparse and Noisy Graphs using Non-Negative Matrix Factorisation."
    % In Proceedings of 22nd ACM International Conference on Information and Knowledge Management, October 2013.       
    %
    % @author: Jeffrey Chan, 2013
    %    
    
    
    
    properties 
        % pre-computed weight matrix.
        m_mWeight;
        m_vWeight;
    end % end of properties

    
    
    methods
        function obj = EucTriFactAdj(mWeight, varargin)
            %
            % constructor
            % mWeight - weight adjustment matrix.
            %
            
            
            obj@EucTriFact(varargin{:});
           
            obj.m_mWeight = mWeight;
            obj.m_vWeight = mWeight(:);
        end
        
        
        function [dis] = distance(obj, mAdj, mImage, mMembership)
            %
            % Compute the approximation distance between mAdj and the
            % approximated blockmodel + regularisation term.
            %

            dis = distance@EucTriFact(obj, mAdj, mImage, mMembership, obj.m_mWeight);
        end % end of distance()
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function [mGrad] = gradient(obj, mAdj, mImage, mMembership)
            %
            % Compute the gradient of approximation function.
            %


            mGrad = 2 * mMembership' * ((mMembership * mImage * mMembership' - mAdj) .* obj.m_mWeight.^2) * mMembership;
   
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
            %

            
            dis = coordDescObj(mAdj, mImage, mMembership, stepSize, basisRow, basisCol, obj.m_mWeight);
           
        end % end of function coordDescDistance
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        
        
        function [mGrad] = posGradient(obj, mAdj, mImage, mMembership)
            %
            % Compute the gradient of the position matrix.
            %
            
            mApproxDiff = (mMembership * mImage * mMembership' - mAdj) .* obj.m_mWeight .* obj.m_mWeight;
            
            mGrad = 2* (mApproxDiff' * mMembership * mImage + mApproxDiff * mMembership * mImage');
        end % end of function posGradient()        
        
    
    end

end % end of class
