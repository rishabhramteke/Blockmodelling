classdef ImageCoordDesc
    %
    % Coordinate descent method for updating the image matrix.
    %
    % If using this code, please kindly cite the following paper:
    % J. Chan, W. Liu, C. Leckie, J. Bailey and K. Ramamohanarao. 
    % "Discovering Latent Blockmodels in Sparse and Noisy Graphs using Non-Negative Matrix Factorisation."
    % In Proceedings of 22nd ACM International Conference on Information and Knowledge Management, October 2013.    
    %
    % @author: Jeffrey Chan, 2014
    %
    
    properties

    end
   
    methods
    
        function obj = ImageCoordDesc(varargin)
         
        end
        
        
        function [mImage] = updateImage(obj, mAdj, mImage, mMembership, fObjective)
            %
            % Updates mImage according to a coordinate descent update rule.
            %
            
             vX0 = fObjective.computeApprox(mAdj, mImage, mMembership);

            % loop through coordinates
            for r = 1 : size(mImage,1)
                for c = 1 : size(mImage,2)
                    vX1 = fObjective.computeX1(mMembership, r, c);
                    [X0dis, X0X1dis, X1dis] = obj.computeTraces(vX0, vX1);
                    % we adjust t until we reach a minimal that satisfies the
                    % constraints on mImage (line search)
                    ftObj = @(t) obj.computeObjectiveFuncOfStepSize(X0dis, X0X1dis, X1dis, t); 
                    [minStepSize, ~] = fminbnd(ftObj, -mImage(r,c), 1 - mImage(r,c));
                    % update mImage
                    mImage(r,c) = minStepSize + mImage(r,c);
                    vX0 = fObjective.updateX0Add(vX0, vX1, minStepSize);
                end % end of inner for

            end % end of outer for           


        end % end of function
    
        
        function [dis] = computeObjectiveFuncOfStepSize(obj, X0dis, X0X1dis, X1dis, stepSize)
                    %
                    % Computes the distance as a function of stepSize, given the other
                    % traces are computed already.
                    %

                    dis = X0dis - 2 * stepSize * X0X1dis + stepSize * stepSize * X1dis;
        end % end of function



        function [ X0dis, X0X1dis, X1dis ] = computeTraces(obj, vX0, vX1)
                    %   
                    % Function that computes the distance of objective according to setting a
                    % particular basis.
                    %

                    % Tr(X0' X0)
                    [~,~,mV0] = find(vX0);
                    X0dis = norm(mV0, 2)^2;

                    % Tr(X0' X1)
                    X0X1dis = full(vX0' * vX1);

                    % Tr(X1' X1)
                    [~,~,mV1] = find(vX1);
                    X1dis = norm(mV1,2)^2;
        end % end of function
        
        
        
    end % end of methods
    
    
end % end of class



