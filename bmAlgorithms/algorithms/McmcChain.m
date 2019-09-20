classdef McmcChain
    %
    % MCMC chain data structure.
    %
    % This assumes two domain and one relation only.
    %
    
    properties
        m_ndomains;
        m_nrelations;
        m_temperature;
        m_chainprob;
        m_mMembership;
    end
    
    methods
        function [obj] = McmcChain(ndomains, nrelations, temperature)
            % 
            % Constructor.
            %
            
            obj.m_ndomains = ndomains;
            obj.m_nrelations = nrelations;
            obj.m_temperature = temperature;
        end
        
        function setMembership(obj, mMembership)
            %
            % Set the membership matrix.
            %
            
            obj.m_mMembership = mMembership;
        end
        
        
    end
    
end % end of class