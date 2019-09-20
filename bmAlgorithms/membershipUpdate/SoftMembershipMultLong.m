classdef SoftMembershipMultLong
    %
    % Long's Multplicative soft membership update function.
    %
    %
    % @author: Jeffrey Chan, 2013
    %

    properties
       m_bNormalisedCol = false; 
       m_bNormalisedRow = false;
       % regularisation parameters
       m_paraL1MembershipReg = 0;
       m_paraL2MembershipReg = 0;       
    end
   
    methods
        
        function obj = SoftMembershipMultLong(varargin)
            if length(varargin) >= 1
               obj.m_bNormalisedCol = varargin{1}; 
            end
            if length(varargin) >= 2
                obj.m_bNormalisedRow = varargin{2};
            end
            if length(varargin) >= 3
                assert(varargin{3} >= 0);
                obj.m_paraL1MembershipReg = varargin{3};
            end
            if length(varargin) == 4
                assert(varargin{4} >= 0);
                obj.m_paraL2MembershipReg = varargin{4};
            end            
        end

            
            
        function [mMembership, mImage] = updateMembership(obj, mAdj, mImage, mMembership)
        %
        % Updates the soft membership, using multiplicatinve update rule.
        %
        %

            mX = mMembership * mImage;
%             mDenom = mX * mMembership' * mX;
%             mNonZero = mDenom ~= 0;
%             mNonim = mAdj * mX;
%             mMembership(mNonZero) = mMembership(mNonZero) .* ((mNonim(mNonZero) ./ mDenom(mNonZero))).^0.25;
%             % this assumes if divide by 0, then membership and mImage is 0 for
%             % that term, so should set mMembership (back) to 0 (0.01 to avoid
%             % inverse singular problems when the image matrix is updated)
%             mMembership(~mNonZero) = 0;
            mMembership = mMembership .* ( (mAdj' * mX) ./ (mX * mMembership' * mX + eps +  obj.m_paraL1MembershipReg + obj.m_paraL2MembershipReg * mMembership) ); 


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
%                 mRowSum = diag(mMembership * ones(size(mMembership,2),1));
%                 mMembership = mRowSum^(-1) * mMembership;
                % do it using for loop, because taking inverse might be singular
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