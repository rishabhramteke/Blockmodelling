classdef ClusterComparisonWrapper
    %
    % Provides a wrapper object around comparison measures.  This wrapper stores
    % a refernce cluster to compare against.
    %
    % @author: Jeffrey Chan, 2014
    %
    
    properties
        m_fComparison = NaN;
        m_mRefMembership = NaN;
    end
    
    
    methods
        function obj = ClusterComparisonWrapper(fComparison, mRefMembership)
            obj.m_fComparison = fComparison;
            obj.m_mRefMembership = mRefMembership;
        end % end of constructor
        
        
        function [val] = compare(obj, mMembership)
            val = obj.m_fComparison(obj.m_mRefMembership, mMembership);
        end % end of compare()
        
    end % end of methods
        
end % end of classdef