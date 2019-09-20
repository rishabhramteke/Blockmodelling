classdef Modularity
    %
    % Computed directed modularity.
    %
    % @author: Jeffrey Chan, 2014
    %
    
    properties
        m_mAdjDeg;
        m_mEdgeNum;
        m_bHarden = true;
    end
    
    
    methods
        function [obj] = Modularity(mAdj, bHarden)
            % compute the adj - deg matrix
            vInDeg = sum(mAdj, 1);
            vOutDeg = sum(mAdj, 2);
            obj.m_mEdgeNum = sum(sum(mAdj));
            obj.m_mAdjDeg = mAdj - (vInDeg' * vOutDeg' / obj.m_mEdgeNum);
            % whether to harden results first before comparison
            obj.m_bHarden = bHarden;
        end
        
        
        function val = modularity(obj, mMembership)
            if obj.m_bHarden
                mNewMembership = discretise(mMembership);
                val = (1 / 2*obj.m_mEdgeNum) * trace(mNewMembership' * obj.m_mAdjDeg * mNewMembership);
            else
                val = (1 / 2*obj.m_mEdgeNum) * trace(mMembership' * obj.m_mAdjDeg * mMembership);
            end
        end
    end % end of methods
    
end % end of class