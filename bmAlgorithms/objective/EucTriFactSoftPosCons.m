classdef EucTriFactSoftPosCons < EucTriFact
    %
    % Eucliean distance with penalty term for soft position constraints to sum
    % to 1.
    %
    %
    % @author: Jeffrey Chan, 2013
    %    
    
    properties
        % position
        m_posConsPara;
    end

    methods
        
        function [obj] = EucTriFactSoftPosCons(varargin)
            obj@EucTriFact();
            
            if ~isempty(varargin)
               assert(length(varargin) == 1);
               obj.m_posConsPara = varargin(1);
            end
        end
        
        
        function [dis] = distance(obj, mAdj, mImage, mMembership, mWeight)
        %
        % Computes the Euclidean distance approximation objective for tri-factorisation.
        %
            dis = distance@EucTriFact(mAdj, mImage, mMebership, mWeight);
            
            % add soft penalty constraint
            dis = dis + obj.m_posConsPara * sum(mMembership * ones(size(mMembership,2),1) - 1);
        end % end of distance()
        
                
                
        function [mGrad] = posGradient(obj, mAdj, mImage, mMembership)
            %
            % Compute the gradient of the position matrix.
            %
            
            mApproxDiff = mMembership * mImage * mMembership' - mAdj;
%             mApproxDiff = mAdj - mMembership * mImage * mMembership';
            
            %mGrad = 2* (mApproxDiff' * mMembership * mImage + mApproxDiff * mMembership * mImage') + obj.m_posConsPara;
            mGrad = obj.posGradientWithApprox(mImage, mMembership, mApproxDiff);
        end % end of function posGradient()
        
        
        function [mGrad] = posGradientWithApprox(obj, mImage, mMembership, mApproxDiff)
            mGrad = 2* (mApproxDiff' * mMembership * mImage + mApproxDiff * mMembership * mImage') + obj.m_posConsPara;
        end
        
    end 

end % end of class
