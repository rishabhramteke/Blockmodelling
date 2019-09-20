classdef BinaryMembershipPenaltyTerm
    
    properties
        m_penatlyWeight = 1;
    end
    
    methods
        function obj = BinaryMembershipPenaltyTerm(varargin)
            if length(varargin) > 1
                obj.m_penatlyWeight = varargin{1};
            end
        end
        
        
        function dis = objective(obj, mMembership, ~, ~)
            %
            % Compute objective.
            %
            
            mApprox = mMembership .* mMembership - mMembership;
            dis = obj.m_penatlyWeight * trace(mApprox' * mApprox);
        end
        
        
        function mGrad = posGradient(obj, ~, ~, mMembership)
           %
           % Compute the position gradient.
           %
           
           mMembershipSq = mMembership .* mMembership;
           mGrad = obj.m_penatlyWeight * ( 4 * mMembershipSq .* mMembership - 6 * mMembershipSq + 2 * mMembership);
        end
    end
    
end % end of classdef