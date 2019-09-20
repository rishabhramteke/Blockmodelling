classdef EucTriFactSoftPosConsAdj < EucTriFactSoftPosCons
    % 
    % Euclidean tri-factorisation with adjusted edge mismatch penalties and soft
    % position constraints.
    %
    %
    % @author: Jeffrey Chan, 2013
    %
    
    properties 
        % pre-computed weight matrix.
        m_mWeight;
    end % end of properties

    
    
    methods
        function obj = EucTriFactSoftPosConsAdj(mWeight, varargin)
            %
            % constructor
            % mWeight - weight adjustment matrix.
            %
            
            
            obj@EucTriFactSoftPosCons(varargin{:});
           
            obj.m_mWeight = mWeight;
        end
        
        
        function [dis] = distance(obj, mAdj, mImage, mMembership)
            %
            % Compute the approximation distance between mAdj and the
            % approximated blockmodel + regularisation term.
            %

            dis = distance@EucTriFactSoftPosCons(obj, mAdj, mImage, mMembership, obj.m_mWeight);
        end % end of distance()
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function [mGrad] = gradient(obj, mAdj, mImage, mMembership)
            %
            % Compute the gradient of approximation function.
            %
            
            % gradient with weight.
            mGrad = 2 * mMembership' * ((mMembership * mImage * mMembership' - mAdj) .* obj.m_mWeight.^2) * mMembership;
        end % end of function gradient()
        
        
        
        
        function [dis] = coordDescDistance(obj, mAdj, mImage, mMembership, stepSize, basisRow, basisCol)
            %
            % Compute the objective value of the coordinate descent expression.
            % TODO: we can reuse most of the calculations again, need to figure
            % out how to store them between calls.
            %
            
            % construct the unit basis
            mUnitBasis = zeros(size(mImage,1), size(mImage,2));
            mUnitBasis(basisRow, basisCol) = 1;
            mUnitBasis = sparse(mUnitBasis);
            
            mX0 = (mAdj - mMembership * mImage * mMembership') .* obj.m_mWeight;
            mX1 = (mMembership * mUnitBasis * mMembership') .* obj.m_mWeight;
            disSq = stepSize * stepSize;
            % Tr(X0' X0)
            dis = trace(mX0' * mX0);
            % Tr(X0' X1)
            dis = dis - 2 * stepSize * trace(mX0' * mX1);
            % Tr(X1' X1)
            dis = dis + disSq * trace(mX1' * mX1);
            
 
        end % end of function coordDescDistance
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        
        
        function [mGrad] = posGradient(obj, mAdj, mImage, mMembership)
            %
            % Compute the gradient of the position matrix.
            %
            
            mApproxDiff = (mMembership * mImage * mMembership' - mAdj) .* obj.m_mWeight .* obj.m_mWeight;
            
            posGradientWithApprox@EucTriFactSoftPosCons(obj, mAdj, mImage, mMembership, mApproxDiff);
        end % end of function posGradient()        
        
    
    end




end % end of class