classdef HardMembershipPenaltyRowConstraint < SoftMembershipMultHardRowConstraint
%
% Ding's Multplicative approach with penatly term for C values that deviate
% from being binary, plus the row constraints to ensure only one entry is
% non-zero.
%
% @author: Jeffrey Chan, 2014
%    

    properties
        
        m_penaltyWeight = 0.01;
        m_penaltyMult = 10;
        m_finalPenaltyWeight = 1;
        m_currIter = 0;
        m_annealingTime = 70;
        m_annealingFact = 0.2;
    end
    
     
    methods
        
        
        function obj = HardMembershipPenaltyRowConstraint(varargin)
        %
        %
        
            obj@SoftMembershipMultHardRowConstraint(varargin{:});        
            
            inParser = inputParser;
            inParser.KeepUnmatched = true;

            % add parameters and set default valuse
            addParameter(inParser, 'initPenaltyWeight', 0.01);
            addParameter(inParser, 'penaltyMult', 10);
            addParameter(inParser, 'finalPenaltyWeight', 1);
            addParameter(inParser, 'initAnnealingTime', 70);
            addParameter(inParser, 'annealingFact', 0.2);
            
            parse(inParser, varargin{:});

            obj.m_penaltyWeight = inParser.Results.initPenaltyWeight;
            obj.m_penaltyMult = inParser.Results.penaltyMult;
            obj.m_finalPenaltyWeight = inParser.Results.finalPenaltyWeight;
            obj.m_annealingTime = inParser.Results.initAnnealingTime;
            obj.m_annealingFact = inParser.Results.annealingFact;


        end % end of function
        
        
        function [mMembership, mImage] = updateMembership(obj, mAdj, mImage, mMembership)
        %
        % Updates the hard membership, using multiplicatinve non-binary penalty terms and row
        % constraints.
        %
        %
            obj.m_currIter = obj.m_currIter + 1;
            
            % see if we need to update the cooling schedule
            if obj.m_currIter >= obj.m_annealingTime
                obj.m_penaltyWeight = obj.m_penaltyWeight * obj.m_penaltyMult;
                obj.m_annealingTime = obj.m_annealingTime + round(obj.m_annealingTime * obj.m_annealingFact);
            end

            posNum = size(mMembership,2);
            mXS = mMembership * mImage;
            mXSt = mMembership * mImage';            
            mMSq = mMembership .* mMembership;
            mMCu = mMembership .* mMSq;            
            
            % if no null-model weights
            if isnan(obj.m_mAdjWeight)
                mMembershipSq = mMembership' * mMembership;

%                                 mPos = mMembership * mImage' * mMembershipSq * mImage + mMembership * mImage * mMembershipSq * mImage' + obj.m_penaltyWeight * (2 * mMCu + mMembership) + obj.m_paraL1MembershipReg + obj.m_paraL2MembershipReg * mMembership;
%                 mNeg = mAdj' * mMembership * mImage + mAdj * mMembership * mImage' + 3 * obj.m_penaltyWeight * mMSq;

                mPos = mXSt * mMembershipSq * mImage + mXS * mMembershipSq * mImage' + obj.m_penaltyWeight * (2 * mMCu + mMembership) + obj.m_paraL1MembershipReg + obj.m_paraL2MembershipReg * mMembership;
                mNeg = mAdj' * mXS + mAdj * mXSt + 3 * obj.m_penaltyWeight * mMSq;

%                 mG = repmat(sum(mMembership./(mPos + eps),2),1, posNum);
%                 mH = repmat(sum(mMembership.*mNeg./(mPos + eps),2),1, posNum);

                mG = repmat(diag(mPos * mMembership'), 1, posNum);
                mH = repmat(diag(mNeg * mMembership'), 1, posNum);

                mMembership = mMembership.*((mNeg + mG)./(mPos + mH + eps)).^0.25;                

%                 mMembership = mMembership.*(mNeg .* mG + 1)./(mPos .* mG + mH + eps);
            % have null-model weights
            else                
%                 mPos = ((mMembership * mImage' * mMembership') .* obj.m_mAdjWeightSq) * mMembership * mImage + ((mMembership * mImage * mMembership') .* obj.m_mAdjWeightSq) * mMembership * mImage' + obj.m_penaltyWeight * (2 * mMCu + mMembership) + obj.m_paraL1MembershipReg + obj.m_paraL2MembershipReg * mMembership;
                mPos = ((mXSt * mMembership') .* obj.m_mAdjWeightSq) * mXS + ((mXS * mMembership') .* obj.m_mAdjWeightSq) * mXSt + obj.m_penaltyWeight * (2 * mMCu + mMembership) + obj.m_paraL1MembershipReg + obj.m_paraL2MembershipReg * mMembership;                                
                mAdjAdj = mAdj .* obj.m_mAdjWeightSq;
                mNeg = mAdjAdj' * mXS + mAdjAdj * mXSt + 3 * obj.m_penaltyWeight * mMSq;

%                 mG = repmat(sum(mMembership./(mPos + eps),2),1, posNum);
%                 mH = repmat(sum(mMembership.*mNeg./(mPos + eps),2),1, posNum);
                
                mG = repmat(diag(mPos * mMembership'), 1, posNum);
                mH = repmat(diag(mNeg * mMembership'), 1, posNum);                
                
                mMembership = mMembership.*((mNeg + mG)./(mPos + mH + eps));   

%                 mMembership = mMembership.*(mNeg .* mG + 1)./(mPos .* mG + mH + eps);                
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