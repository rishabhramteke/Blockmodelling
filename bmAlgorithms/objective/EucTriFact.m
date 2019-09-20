classdef EucTriFact
    %
    % @author: Jeffrey Chan, 2013
    %    
    
    properties
       m_fComputeCoordDesc = ComputeCoordDescObjective; 
       % Other extra terms
       m_fComputeExtraTerms = NaN;
    end
    

    methods
        
        function obj = EucTriFact(varargin)
            inParser = inputParser;
            inParser.KeepUnmatched = true;
            
            addParameter(inParser, 'computeExtraTerms', NaN);

            parse(inParser, varargin{:});

            obj.m_fComputeExtraTerms = inParser.Results.computeExtraTerms;            
        end
        
        function [dis] = distance(obj, mAdj, mImage, mMembership, varargin)
        %
        % Computes the Euclidean distance approximation objective for tri-factorisation.
        %

            %objVal = norm(mAdj - mMembership * mImage * mMembership',2);
%             mApprox = mMembership * mImage * mMembership';
%             mApprox = mAdj - mApprox;
            mApprox = mAdj - mMembership * mImage * mMembership';
            
            % see if there is a weight matrix
            if length(varargin) == 1
                mWeight = varargin{1};
                mApprox = mApprox .* mWeight;
            end
            
            [~,~,mV] = find(mApprox); 
            dis = norm(mV, 2)^2;
%             dis = trace(mApprox * mApprox');

            % convert to number if in sparse representation.
            if issparse(dis)
                dis = full(dis);
            end
            
            if ~isnan(obj.m_fComputeExtraTerms)
                dis = dis + obj.m_fComputeExtraTerms.objective(mMembership, mImage, mAdj);
            end
        end % end of distance()
        
        
        function [mGrad] = gradient(~, mAdj, mImage, mMembership)
            %
            % Compute the gradient of the image matrix.
            %
%             
            mMemTerm = mMembership' * mMembership;
            
%             mApprox = mMembership * mImage * mMembership' - mAdj;
%             mGrad = 2 * mMembership' * mApprox * mMembership;
                        
            % add all the pieces intogether
            mGrad =  2 * (mMemTerm * mImage * mMemTerm - mMembership' * mAdj * mMembership);
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
        
        
        
        
        function [mGrad] = posGradient(obj, mAdj, mImage, mMembership)
            %
            % Compute the gradient of the position matrix.
            %
            
%             mApproxDiff = mMembership * mImage * mMembership' - mAdj;
% %             mApproxDiff = mAdj - mMembership * mImage * mMembership';
%             
%             mGrad = 2* (mApproxDiff' * mMembership * mImage + mApproxDiff * mMembership * mImage');
                mSame = mMembership' * mMembership;
                mGrad = 2 * (mMembership * mImage' * mSame * mImage + mMembership * mImage * mSame * mImage' -...
                    mAdj' * mMembership * mImage + mAdj * mMembership * mImage');
                
                if ~isnan(obj.m_fComputeExtraTerms)
                    mGrad = mGrad + obj.posGradientExtraTerm(mAdj, mImage, mMembership);
                end
        end % end of function posGradient()
        
        
        function [mGrad] = posGradientExtraTerm(obj, mAdj, mImage, mMembership)
            % 
            % Compute the gradient for the extra terms.
            %
            
            mGrad = NaN;
            if ~isnan(obj.m_fComputeExtraTerms)
                mGrad = obj.m_fComputeCoordDesc.posGradient(mAdj, mImage, mMembership);
            end
        end
        
    end 

end % end of class
