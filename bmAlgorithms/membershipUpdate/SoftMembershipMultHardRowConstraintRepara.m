classdef SoftMembershipMultHardRowConstraintRepara
    %
    % Multplicative soft membership update function with hard row constraints, 
    % using the reparameterisation technique.
    %
    % Introduce row normalisation of vertex membership as well as L1 and L2
    % regularisation as options.
    %
    %
    % @author: Jeffrey Chan, 2014
    %    
    
    properties
       m_bNormalisedCol = false; 
       % regularisation parameters
       m_paraL1MembershipReg = 0;
       m_paraL2MembershipReg = 0;
    end
    
    
    methods
        
        
        function obj = SoftMembershipMultHardRowConstraintRepara(varargin)
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
                assert(varargin{2} >= 0);
                obj.m_paraL1MembershipReg = varargin{2};
            end
            if length(varargin) ==3
                assert(varargin{3} >= 0);
                obj.m_paraL2MembershipReg = varargin{3};
            end
           
        end % end of function
        
        
        function [mMembership, mImage] = updateMembership(obj, mAdj, mImage, mMembership)
        %
        % Updates the soft membership, using multiplicatinve update rule.
        %
        %

            posNum = size(mMembership,2);
            
            mMembershipSq = mMembership' * mMembership;
            mPos = mMembership * mImage' * mMembershipSq * mImage + mMembership * mImage * mMembershipSq * mImage' + obj.m_paraL1MembershipReg + obj.m_paraL2MembershipReg * mMembership;
            mNeg = mAdj' * mMembership * mImage + mAdj * mMembership * mImage';
        
            mG = repmat(diag(mPos * mMembership'), 1, posNum);
            mH = repmat(diag(mNeg * mMembership'), 1, posNum);
        
            mMembership = mMembership.*(mNeg + mG)./(mPos + mH + eps);
        
            % normalisation
            mMembership = mMembership ./ repmat(sum(mMembership,2), 1, posNum);
            
            mMembership = max(eps, mMembership);
        
            % normalise each column to unit length
            if obj.m_bNormalisedCol
                % column sum
                mD = diag(sum(mMembership, 1));
                mMembership = mMembership / mD;
                mImage = mD * mImage * mD;
%                 mMembership = normc(mMembership);
            end
            
        end % end of function
        
        
    end % end of methods
    
end