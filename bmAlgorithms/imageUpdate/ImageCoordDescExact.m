classdef ImageCoordDescExact
    %
    % @author: Jeffrey Chan, 2014
    %
    
    properties
        % set to true if small amount of gaussian noise to be added if B set to
        % 0 or 1
        m_bAddNoise = false;
        
        m_l1ImageReg = 0;
    end
   
    methods
    
        function obj = ImageCoordDescExact(varargin)
            if length(varargin) >= 1
                obj.m_bAddNoise = varargin{1};
            end
            if length(varargin) == 2
                obj.m_l1ImageReg = varargin{2};
            end
        end
        
        
        function [mImage] = updateImage(obj, mAdj, mImage, mMembership, fDistanceFunc)
            %
            % Updates mImage according to a coordinate descent update rule, with exact step
            % size calculated (no need for line search).
            % Should only be used for the euclidean distance.
            %
            % 
            %
            
%                 mUnitBasis = zeros(size(mImage,1), size(mImage,2));

                % euclidean
                mX0 = mAdj - mMembership * mImage * mMembership';
                vX0Trans = mX0(:)';

                % loop through coordinates
                for r = 1 : size(mImage,1)
                    for c = 1 : size(mImage,2)
%                         mUnitBasis(r,c) = 1;
%                         mX1a = mMembership * mUnitBasis * mMembership';
                        % outer product
                        mX1 = mMembership(:,r) * mMembership(:,c)';
                        vX1 = mX1(:);
                        % computing trace(mX1' * mX1)
                        denom = full(vX1' * vX1);
%                         denom2 = trace(mX1' * mX1);
%                         denom3 = trace(mX1a' * mX1a);
            %             denom = full(vX1' * vX1);
            %             fprintf('denom = %.2f\n', denom);
            %             fprintf('trace() = %.2f\n', trace(mX1' * mX1));

                        nomin = full(vX0Trans * vX1) - 0.5 * obj.m_l1ImageReg;
%                         nomin2 = trace(mX0'*mX1);
            %             fprintf('nomin = %.2f\n', nomin);
            %             fprintf('trace() = %.2f\n', trace(mX0' * mX1));            
                        % compute optimal step size
            %             altStepSize = trace(mX0' * mX1)/trace(mX1' * mX1);
                        if denom == 0
                            denom = eps;
                        end
                        
                        altStepSize = nomin / denom;
                        
                        
                        if altStepSize > 0
                            minStepSize = min(1 - mImage(r,c), altStepSize);
                        else
                            minStepSize = max(-mImage(r,c), altStepSize);
                        end
                        % update mImage
                        mImage(r,c) = minStepSize + mImage(r,c);
                        
                        if obj.m_bAddNoise
                            if mImage(r,c) == 0
                                if rand > 0.5
                                    mImage(r,c) = max(0, normrnd(0.0001, 0.00003));
                                end
                            elseif mImage(r,c) == 1
                                if rand > 0.5
                                    mImage(r,c) = min(1, normrnd(0.9999, 0.00003));
                                end
                            end
                        end 
%                         mUnitBasis(r,c) = 0;
                    end % end of inner for

                end % end of outer for       

        end % end of function
        
    end % end of methods
    
    
end % end of class



