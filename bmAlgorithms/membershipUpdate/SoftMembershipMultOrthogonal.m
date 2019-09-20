classdef SoftMembershipMultOrthogonal
    %
    % Multplicative soft membership update function, orthogonal constraint.
    % From "Quadratic Non-negative matrix factorisation paper" of Yang and Oja.
    %
    %
    % @author: Jeffrey Chan, 2014
    %    
    
    properties

    end
    
    
    methods
        
        
        function obj = SoftMembershipMultOrthogonal(varargin)

        end % end of function
        
        
        function [mMembership, mImage] = updateMembership(~, mAdj, mImage, mMembership)
        %
        % Updates the soft membership, using multiplicatinve update rule.
        %
        %
%             fprintf('%s\n', 'membership');
%             if ~isreal(mMembership) 
%                 display(mMembership);
%             end
%             
%             if ~isreal(mImage) 
%                 display(mImage);
%             end        
        
            mMemSqInv = mMembership' * mMembership;
            mPos = mMembership * mImage' * mMemSqInv * mImage + mMembership * mImage * mMemSqInv * mImage';
            mNeg = mAdj' * mMembership * mImage + mAdj * mMembership * mImage';            
            
            mMemSq = mMembership * mMembership';
            mMembership = mMembership .* ( ( (mNeg + mMemSq * mPos) ./ (mPos + mMemSq * mNeg + eps) ).^0.5);
            
        end % end of function
        
        
    end % end of methods
    
end