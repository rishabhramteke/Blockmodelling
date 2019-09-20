classdef ImageMultDing
    %
    % Ding's Multplicative soft image update function.
    %
    % Introduce row normalisation of vertex membership as well as L1 and L2
    % regularisation as options.
    %
    %
    % @author: Jeffrey Chan, 2014
    %    
    
    properties
       % regularisation parameters
       m_paraL1MembershipReg = 0;
       m_paraL2MembershipReg = 0;
    end
    
    
    methods
        
        
        function obj = ImageMultDing(varargin)
        %
        % INPUT:
        % varargin{1}:
        % varargin{2}:
        %
              
            if length(varargin) >= 1
                assert(varargin{1} >= 0);
                obj.m_paraL1MembershipReg = varargin{1};
            end
            if length(varargin) == 2
                assert(varargin{2} >= 0);
                obj.m_paraL2MembershipReg = varargin{2};
            end
        end % end of function
        
        
        function [mImage] = updateImage(obj, mAdj, mImage, mMembership, fDistanceFunc)
        %
        % Using multiplicatinve update rule.
        %
        %
        
%             fprintf('%s\n', 'image');
%             if ~isreal(mMembership) 
%                 display(mMembership);
%             end
%             
%             if ~isreal(mImage) 
%                 display(mImage);
%             end             

            mXtX = mMembership' * mMembership;
%       mDenom = (mXtX * mImage * mXtX);
%       mNomin = (mMembership' * mAdj * mMembership);
%       % avoid divide by zero errors
%       mNonZero = mDenom ~= 0;
%       mImage(mNonZero) = mImage(mNonZero) .* (mNomin(mNonZero) ./ mDenom(mNonZero));
%       mImage(~mNonZero) = 0;
    
%       mFactor = (mMembership' * mAdj * mMembership) ./ (mXtX * mImage * mXtX);
%       % remove any divide by zero errors
%       mFactor(isnan(mFactor)) = 0;
%       mFactor(isinf(mFactor)) = 0;



            % update image
            mImage = mImage .* (mMembership' * mAdj * mMembership) ./ (mXtX * mImage * mXtX + eps + obj.m_paraL1MembershipReg + obj.m_paraL2MembershipReg * mImage);
  
            % help algorithm not to get stuck on 0 entries
            mImage = max(eps, mImage);            
        end % end of function
        
        
    end % end of methods
    
end