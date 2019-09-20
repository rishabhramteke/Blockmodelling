classdef SoftMembershipMultLongIDiv
    %
    % Long's Multplicative soft membership update function for I-divergence
    %
    %
    % @author: Jeffrey Chan, 2014
    %

    properties
       m_bNormalisedCol = false; 
       m_bNormalisedRow = false;
    end
   
    methods
        
        function obj = SoftMembershipMultLong(varargin)
            if length(varargin) >= 1
               obj.m_bNormalisedCol = varargin{1}; 
            end
            if length(varargin) == 2
                obj.m_bNormalisedRow = varargin{2};
            end
        end

            
            
        function [mMembership, mImage] = updateMembership(obj, mAdj, mImage, mMembership)
        %
        % Updates the soft membership, using multiplicatinve update rule.
        %
        %

            mMembership = mMembership .* ((mAdj ./ (mMembership * mImage * mMembership')) * mMembership * mImage) ./ (ones(size(mAdj,1), size(mAdj,1)) * mMembership * mImage)).^(0.5):


            % normalise each column to unit length
            if obj.m_bNormalisedCol
                % column sum
                vDiagsum = sum(mMembership, 1);
                mD = diag(vDiagsum);
                mDInv = diag(1 ./ vDiagsum);
                mMembership = mMembership * mDInv;
                mImage = mD * mImage * mD;
%                 mMembership = normc(mMembership);
            end
            
            
            if obj.m_bNormalisedRow
                mRowSum = diag(mMembership * ones(size(mMembership,2),1));
                mMembership = mRowSum^(-1) * mMembership;
            end
        end % end of function
        
        
    end % end of methods
    
end