classdef SoftMembershipMultDing
    %
    % Ding's Multplicative soft membership update function.
    %
    % Introduce row normalisation of vertex membership as well as L1 and L2
    % regularisation as options.
    %
    %
    % @author: Jeffrey Chan, 2013
    %    
    
    properties
       m_bNormalisedCol = false; 
       m_bAdhocNormalisedRow = false;
       % regularisation parameters
       m_paraL1MembershipReg = 0;
       m_paraL2MembershipReg = 0;
    end
    
    
    methods
        
        
        function obj = SoftMembershipMultDing(varargin)
        %
        % INPUT:
        % varargin{1}:
        % varargin{2}:
        % varargin{3}:
        %
        
            if length(varargin) >= 1
               obj.m_bNormalisedCol = varargin{1}; 
            end
            if length(varargin) >= 2
                obj.m_bAdhocNormalisedRow = varargin{2};
            end            
            if length(varargin) >= 3
                assert(varargin{3} >= 0);
                obj.m_paraL1MembershipReg = varargin{3};
            end
            if length(varargin) == 4
                assert(varargin{4} >= 0);
                obj.m_paraL2MembershipReg = varargin{4};
            end
            
        end % end of function
        
        
        function [mMembership, mImage] = updateMembership(obj, mAdj, mImage, mMembership)
        %
        % Updates the soft membership, using multiplicatinve update rule.
        %
        %

            mXS = mMembership * mImage;
            mXSt = mMembership * mImage';
%             mDenom = (mXS * (mMembership') * mXSt + mXSt * (mMembership') * mXS);
%             mNomin = (mAdj' * mXS + mAdj * mXSt);
%             mNonZero = mDenom ~= 0;
%             mMembership(mNonZero) = mMembership(mNonZero) .* (mNomin(mNonZero) ./ mDenom(mNonZero)).^0.25;
%             mMembership(~mNonZero) = 0;
            mMembership = mMembership .* ( (mAdj' * mXS + mAdj * mXSt) ./ (mXS * mMembership' * mXSt + mXSt * mMembership' * mXS + eps + obj.m_paraL1MembershipReg + obj.m_paraL2MembershipReg * mMembership) ).^0.25; 
            
            
            % normalise each column to unit length
            if obj.m_bNormalisedCol
                % column sum
                mD = diag(sum(mMembership, 1));
                mDInv = inv(mD);
                mMembership = mMembership * mDInv;
                mImage = mD * mImage * mD;
%                 mMembership = normc(mMembership);
            end
            
            if obj.m_bAdhocNormalisedRow
%                 mRowSum = diag(mMembership * ones(size(mMembership,2),1));
%                 mMembership = mRowSum^(-1) * mMembership;
                vRowSum = sum(mMembership, 2);
                for r = 1 : length(vRowSum)
                    if vRowSum(r) > eps
                        mMembership(r,:) = mMembership(r,:) / vRowSum(r);
                    end
                end                
            end            
        end % end of function
        
        
    end % end of methods
    
end