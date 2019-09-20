classdef IDivTriFact
    %
    % I-Divergence for tri-factorisation formulation.
    %
    % @author: Jeffrey Chan, 2014
    %    
    
    properties
       m_fComputeCoordDesc = ComputeCoordDescObjective; 
    end
    

    methods
        function [dis] = distance(~, mAdj, mImage, mMembership)
        %
        % Computes the I-Divergence approximation objective for tri-factorisation.
        %

            mApprox = mMembership * mImage * mMembership';
            
            dis = mApprox .* log((mApprox + eps) ./ (mAdj + eps)) - mApprox + mAdj;
            
            % convert to number if in sparse representation.
            if issparse(dis)
                dis = full(dis);
            end
        end % end of distance()
        
        
        function [mGrad] = gradient(~, mAdj, mImage, mMembership)
            %
            % Compute the gradient of the image matrix.
            %
%             
%             mFirst = -2 * mMembership' * mAdj * mMembership;
%             mMemTerm = mMembership' * mMembership;
%             mSecond = 2 * (mMemTerm * mImage * mMemTerm);
            
            mApprox = mMembership * mImage * mMembership' - mAdj;
            mGrad = 2 * mMembership' * mApprox * mMembership;
                        
            % add all the pieces intogether
%             mGrad = mFirst + mSecond;
        end % end of function gradient()
        
        function vX0 = computeApprox(obj, mAdj, mImage, mMembership)
            vX0 = obj.m_fComputeCoordDesc.computeApprox(mAdj, mImage, mMembership);
        end % end of function
        
        

        function [ vX1 ] = computeX1(obj, mMembership, basisRow, basisCol)
           vX1 =  obj.m_fComputeCoordDesc.computeX1(mMembership, basisRow, basisCol);
        end
        
        
        function vX0 = updateX0Add(obj, vX0, mX1, stepSize)
            vX0 = obj.m_fComputeCoordDesc.updateX0Add(vX0, mX1, stepSize);
        end
        
        
        function [dis, vX1] = coordDescDistance(obj, mAdj, mImage, mMembership, stepSize, basisRow, basisCol)
            %
            % Compute the objective value of the coordinate descent expression
            % for the image matrix.
            % TODO: we can reuse most of the calculations again, need to figure
            % out how to store them between calls.
            %
            
            % construct the unit basis
%             mUnitBasis = zeros(size(mImage,1), size(mImage,2));
%             mUnitBasis(basisRow, basisCol) = 1;
%             mUnitBasis = sparse(mUnitBasis);
            
            % Tr(X0' X0)
%             dis = obj.distance(mAdj, mImage, mMembership);
%             mX1 = mMembership * mUnitBasis * mMembership';
%             dis = dis - 2 * stepSize * trace((mAdj - mMembership * mImage * mMembership')' * mX1);
%             disSq = stepSize * stepSize;
%             dis = dis + disSq * trace(mX1' * mX1);
            
            
            [dis, vX1] = obj.m_fComputeCoordDesc.computeObjective(mAdj, mImage, mMembership, stepSize, basisRow, basisCol);
        end % end of function coordDescDistance      
        
        
        
        
        function [mGrad] = posGradient(~, mAdj, mImage, mMembership)
            %
            % Compute the gradient of the position matrix.
            %
            
            mApproxDiff = mMembership * mImage * mMembership' - mAdj;
%             mApproxDiff = mAdj - mMembership * mImage * mMembership';
            
            mGrad = 2* (mApproxDiff' * mMembership * mImage + mApproxDiff * mMembership * mImage');
        end % end of function posGradient()
        
    end 

end % end of class
