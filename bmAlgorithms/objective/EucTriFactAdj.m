classdef EucTriFactAdj < EucTriFact
    % 
    % Euclidean tri-factorisation and has adjusted
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
            % TODO: we can reuse most of the calculations again, need to figure
            % out how to store them between calls.
            %

            
            % construct the unit basis
%             mUnitBasis = zeros(size(mImage,1), size(mImage,2));
%             mUnitBasis(basisRow, basisCol) = 1;
%             mUnitBasis = sparse(mUnitBasis);
            
%             mX0 = (mAdj - mMembership * mImage * mMembership') .* obj.m_mWeight;
% %             mX1 = (mMembership * mUnitBasis * mMembership') .* obj.m_mWeight;
%             mX1 = (mMembership(:,basisRow) * mMembership(:,basisCol)') .* obj.m_mWeight;
% %             vX0 = mX0(:);
% %             vX1 = mX1(:);
%             disSq = stepSize * stepSize;
%             % Tr(X0' X0)
%             [~,~,mV0] = find(mX0); 
%             dis = norm(mV0, 2)^2;
% %             dis = full(vX0' * vX0);
% %             dis = trace(mX0' * mX0);
%             % Tr(X0' X1)
%             dis = dis - 2 * stepSize * trace(mX0' * mX1);
% %             dis = dis - 2 * stepSize * full(vX0' * vX1);
%             % Tr(X1' X1)
% %             dis = dis + disSq * trace(mX1' * mX1);
%             % trace(mX1' * mX1) = sum(mMembership(*,r) outer_product
%             % mMembership(*.c)).^2)
%             [~,~,mV1] = find(mX1);
%             
% %             dis = dis + disSq * full(vX1' * vX1);
%             dis = dis + disSq * norm(mV1,2)^2;
            dis = coordDescObj(mAdj, mImage, mMembership, stepSize, basisRow, basisCol, obj.m_mWeight);
           
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
