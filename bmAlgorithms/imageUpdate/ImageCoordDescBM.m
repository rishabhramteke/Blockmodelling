classdef ImageCoordDescBM
    %
    % @author: Jeffrey Chan, 2014
    %
    
    properties

    end
   
    methods
    
        function obj = ImageCoordDescBM()
           
        end
        
        
        function [mImage] = updateImage(obj, mAdj, mImage, mMembership, fDistanceFunc)
            %
            % Updates mImage according to a coordinate descent update rule.
            %
            
            vX0 = fDistanceFunc.computeApprox(mAdj, mImage, mMembership);

            % loop through coordinates
            for r = 1 : size(mImage,1)
                for c = 1 : size(mImage,2)
                    vX1 = fDistanceFunc.computeX1(mMembership, r, c);
                    [X0dis, X0X1dis, X1dis] = obj.computeTraces(vX0, vX1);
                    % we adjust t until we reach a minimal that satisfies the
                    % constraints on mImage (line search)
        %             ftObj = @(t) fObjective.coordDescDistance(mAdj, mImage, mMembership, t, r, c); 
                    ftObj = @(t) obj.computeObjectiveFuncOfStepSize(X0dis, X0X1dis, X1dis, r, c, mAdj, mImage, fDistanceFunc, t); 
%                     ftObj = @(t) obj.computeObjectiveFuncOfStepSize(X0dis, X0X1dis, X1dis, t); 
                    [minStepSize, ~] = fminbnd(ftObj, -mImage(r,c), 1 - mImage(r,c));
                    % update mImage
                    mImage(r,c) = minStepSize + mImage(r,c);
                    vX0 = fDistanceFunc.updateX0Add(vX0, vX1, minStepSize);
                end % end of inner for

            end % end of outer for       


        end % end of function
    
        
        function [dis] = computeObjectiveFuncOfStepSize(~, X0dis, X0X1dis, X1dis, basisRow, basisCol, mAdj, mImage, fObjective, stepSize)
                    %
                    % Computes the distance as a function of stepSize, given the other
                    % traces are computed already.
                    %

                    dis = X0dis - 2 * stepSize * X0X1dis + stepSize * stepSize * X1dis +...
                        fObjective.computeImageRegDis(mAdj, mImage, stepSize, basisRow, basisCol);
        end % end of function



        function [ X0dis, X0X1dis, X1dis ] = computeTraces(~, vX0, vX1)
                    %   
                    % Function that computes the distance of objective according to setting a
                    % particular basis.
                    %
                    % mAdj -    Adjacency matrix.
                    % mImage -  Image matrix.
                    % mMembership - Vertex position membership matrix.
                    % stepSize - Size to increase basis by.
                    % basisRow
                    % basisCol
                    % varargin - If set, should be the adjustment weighting matrix.
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



