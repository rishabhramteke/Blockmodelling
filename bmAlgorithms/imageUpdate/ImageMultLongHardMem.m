classdef ImageMultLongHardMem
    %
    % Long's Multplicative image update function, for hard membership.
    %
    % Updates mImage according to KKT condition derived multiplicative update
    % rule.
    %
    % Bo Long formulation for Hard cluster membership
    %
    %
    % @author: Jeffrey Chan, 2014
    %   
    
    properties

    end % end of properties
    
    
    methods
        
        function obj = ImageMultLongHardMem()

        end % end of function
        
        
        function [mImage] = updateImage(~, mAdj, mImage, mMembership, fDistanceFunc)
        %
        % Using multiplicatinve update rule.
        %
        %
        
            % mMembership
            % mMembership' * mMembership
            % Using pseudo inverse, to prevent divide by zero issues
            % mSqInv = sparse(inv(mMembership' * mMembership));
            mSqInv = sparse(pinv(full(mMembership' * mMembership)));

            % shouldn't need this if using pseudo-inverse

            % if we get Nan values for mSqInv, it means some row and columsn of mMembership
            % are all zeros or close to all zeros, which could occur due to computation
            % accuracy limitations.  We set these values to 0 instead.
            % mSqInv(isnan(mSqInv)) = 0;

            mImage = mSqInv * mMembership' * mAdj * mMembership * mSqInv;
  
        end % end of function
        
        
    end % end of methods
    
end