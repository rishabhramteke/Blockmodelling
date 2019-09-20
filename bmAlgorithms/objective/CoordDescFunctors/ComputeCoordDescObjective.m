classdef ComputeCoordDescObjective
    %
    % @author: Jeffrey Chan, 2013
    %    

    properties
%         m_vX0 = NaN;
%         m_mX1 = NaN;
    end % end of properties
    
    methods
    
        function obj = ComputeCoordDescObjective()
            % Constructor
    
        end
        
        
        function vX0 = computeApprox(obj, mAdj, mImage, mMembership)
            %
            % Computes the mX0 matrix approximate to A and stores it internally.
            %
                

                mX0 = mAdj - mMembership * mImage * mMembership';
                vX0 = mX0(:);

        end % end of function
        
        
        function vX0 = computeWeightedApprox(obj, mAdj, mImage, mMembership, mWeight)
            %
            % Computes the mX0 matrix approximate to A and stores it internally.
            %
            

                mX0 = (mAdj - mMembership * mImage * mMembership') .* mWeight;
                vX0 = mX0(:);
        end % end of function        
    
    

        function vX0 = updateX0Add(obj, vX0, vX1, stepSize)
            vX0 = vX0 - stepSize * vX1;
        end
        
        
        function vX0 = updateWeightedX0Add(obj,  vX0, vX1, stepSize, vWeight)
            vX0 = vX0 - stepSize * (vX1 .* vWeight);
        end
        
            
        function [ vX1 ] = computeX1(obj, mMembership, basisRow, basisCol)
            mX1 = mMembership(:,basisRow) * mMembership(:,basisCol)';
            vX1 = mX1(:);
        end
        
        function [ vX1 ] = computeWeightedX1(obj, mMembership, basisRow, basisCol, mWeight)
            mX1 = (mMembership(:,basisRow) * mMembership(:,basisCol)') .* mWeight;
            vX1 = mX1(:);
        end
        
%     
%                             mWeight = varargin{1};
% %                     mX0 = (mAdj - mMembership * mImage * mMembership') .* mWeight;
% %             mX1 = (mMembership * mUnitBasis * mMembership') .* obj.m_mWeight;
%                     mX1 = (mMembership(:,basisRow) * mMembership(:,basisCol)') .* mWeight;
        

        
        

    end % end of methods

end % end of class