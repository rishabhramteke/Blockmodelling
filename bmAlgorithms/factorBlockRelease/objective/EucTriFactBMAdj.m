classdef EucTriFactBMAdj < EucTriFactBM
    % 
    % Euclidean tri-factorisation frobenius objective with regularisation term and has adjusted
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
            
            % this is the term that is different from EucTriFactBM
            mTerm = 2 * mMembership' * ((mMembership * mImage * mMembership' - mAdj) .* obj.m_mWeight.^2) * mMembership;

            
            mRegGrad = obj.imageRegGradient(mImage);
               
            % add all the pieces intogether
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
            %
            
                    
            dis = coordDescObj(mAdj, mImage, mMembership, stepSize, basisRow, basisCol, obj.m_mWeight);            
                   
            dis = dis + obj.m_regPara *  obj.imageRegCoordDesc(mImage, stepSize, basisRow, basisCol);            
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