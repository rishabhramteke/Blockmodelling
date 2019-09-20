classdef SoftMembershipMultHardRowConstraint
    %
    % Multplicative soft membership update function with hard row constraints.
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
       m_mAdjWeight;
       m_mAdjWeightSq = NaN;
    end
    
    
    methods
        
        
        function obj = SoftMembershipMultHardRowConstraint(varargin)
        %
        % INPUT:
        % varargin{1}:
        % varargin{2}:
        % varargin{3}:
        %
            inParser = inputParser;
            inParser.KeepUnmatched = true;

            % add parameters and set default valuse
            addParameter(inParser, 'normaliseCol', false);
            addParameter(inParser, 'l1RegWeight', 0);
            addParameter(inParser, 'l2RegWeight', 0);
            addParameter(inParser, 'adjWeight', NaN);
            

            parse(inParser, varargin{:});

            obj.m_bNormalisedCol = inParser.Results.normaliseCol;
            obj.m_paraL1MembershipReg = inParser.Results.l1RegWeight;
            obj.m_paraL2MembershipReg = inParser.Results.l2RegWeight;
            obj.m_mAdjWeight = inParser.Results.adjWeight;
            if ~isnan(obj.m_mAdjWeight)
                obj.m_mAdjWeightSq = obj.m_mAdjWeight .* obj.m_mAdjWeight;
            end
                   
        end % end of function
        
        
        function [mMembership, mImage] = updateMembership(obj, mAdj, mImage, mMembership)
        %
        % Updates the soft membership, using multiplicatinve update rule.
        %
        %

            posNum = size(mMembership,2);
            
            % copied form Oscar's code
            if isnan(obj.m_mAdjWeight)
                mMembershipSq = mMembership' * mMembership;
                GRAD_POS = mMembership * mImage' * mMembershipSq * mImage + mMembership * mImage * mMembershipSq * mImage' + obj.m_paraL1MembershipReg + obj.m_paraL2MembershipReg * mMembership;
                GRAD_NEG = mAdj' * mMembership * mImage + mAdj * mMembership * mImage';

                mG = repmat(sum(mMembership./(GRAD_POS + eps),2),1, posNum);
                mH = repmat(sum(mMembership.*GRAD_NEG./(GRAD_POS + eps),2),1, posNum);

                mMembership = mMembership.*(GRAD_NEG .* mG + 1)./(GRAD_POS .* mG + mH + eps);
            else                
%                 ((mMembership * mImage' * mMembership') .* obj.m_mAdjWeightSq) * mMembership * mImage;
%                 ((mMembership * mImage * mMembership') .* obj.m_mAdjWeightSq) * mMembership * mImage';
%                 obj.m_paraL1MembershipReg + obj.m_paraL2MembershipReg * mMembership;
                mPos = ((mMembership * mImage' * mMembership') .* obj.m_mAdjWeightSq) * mMembership * mImage + ((mMembership * mImage * mMembership') .* obj.m_mAdjWeightSq) * mMembership * mImage' + obj.m_paraL1MembershipReg + obj.m_paraL2MembershipReg * mMembership;
                mAdjAdj = mAdj .* obj.m_mAdjWeightSq;
                mNeg = mAdjAdj' * mMembership * mImage + mAdjAdj * mMembership * mImage';

                mG = repmat(sum(mMembership./(mPos + eps),2),1, posNum);
                mH = repmat(sum(mMembership.*mNeg./(mPos + eps),2),1, posNum);

                mMembership = mMembership.*(mNeg .* mG + 1)./(mPos .* mG + mH + eps);                
            end
            
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