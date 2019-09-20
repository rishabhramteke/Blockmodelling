classdef EucTriFact
    %
    % Euclidean tri-factorisation frobenius objective.
    %
    % If using this code, please kindly cite the following paper:
    % J. Chan, W. Liu, C. Leckie, J. Bailey and K. Ramamohanarao. 
    % "Discovering Latent Blockmodels in Sparse and Noisy Graphs using Non-Negative Matrix Factorisation."
    % In Proceedings of 22nd ACM International Conference on Information and Knowledge Management, October 2013.       
    %
    % @author: Jeffrey Chan, 2013
    %    
    
    properties
       m_fComputeCoordDesc = ComputeCoordDescObjective; 
    end
    

    methods
        function [dis] = distance(~, mAdj, mImage, mMembership, mWeight)
        %
        % Computes the Euclidean distance approximation objective for tri-factorisation.
        %


            mApprox = mAdj - mMembership * mImage * mMembership';
            
            % see if there is a weight matrix
            if nargin == 5
                mApprox = mApprox .* mWeight;
            end
            
            [~,~,mV] = find(mApprox); 
            dis = norm(mV, 2)^2;

            % convert to number if in sparse representation.
            if issparse(dis)
                dis = full(dis);
            end
        end % end of distance()
        
        
        function [mGrad] = gradient(~, mAdj, mImage, mMembership)
            %
            % Compute the gradient of the image matrix.
            %
           
            mMemTerm = mMembership' * mMembership;
            
                
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
            %
            
            
            [dis, vX1] = obj.m_fComputeCoordDesc.computeObjective(mAdj, mImage, mMembership, stepSize, basisRow, basisCol);
        end % end of function coordDescDistance      
        
        
        
        
        function [mGrad] = posGradient(~, mAdj, mImage, mMembership)
            %
            % Compute the gradient of the position matrix.
            %
            
                mSame = mMembership' * mMembership;
                mGrad = 2 * (mMembership * mImage' * mSame * mImage + mMembership * mImage * mSame * mImage' -...
                    mAdj' * mMembership * mImage + mAdj * mMembership * mImage');
        end % end of function posGradient()
        
    end 

end % end of class
