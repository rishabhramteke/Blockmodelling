classdef ImageCoordDescBMExact
    %
    % Update the image using coordinate descent and with the cubic BM function. 
    %
    % @author: Jeffrey Chan, 2014
    %
    
    properties
        m_l1ImageReg = 0;
    end
   
    methods
    
        function obj = ImageCoordDescBMExact(varargin)
            if length(varargin) == 1
                obj.m_l1ImageReg = varargin{1};
            end
         
        end
        
        
        function [mImage] = updateImage(obj, mAdj, mImage, mMembership, fDistanceFunc)
            %
            % Updates mImage according to a coordinate descent update rule.
            %
            
            % TODO incrementally update the coordinate components.
            
            vX0 = fDistanceFunc.computeApprox(mAdj, mImage, mMembership);

            % loop through coordinates
            for r = 1 : size(mImage,1)
                for c = 1 : size(mImage,2)
                    vX1 = fDistanceFunc.computeX1(mMembership, r, c);
                    [X0dis, X0X1dis, X1dis] = obj.computeTraces(vX0, vX1);
                    % we adjust t until we reach a minimal that satisfies the
                    % constraints on mImage (line search)
                    
                    
        %             ftObj = @(t) fObjective.coordDescDistance(mAdj, mImage, mMembership, t, r, c); 
                    minStepSize = obj.computeOptimalStepSize(X0dis, X0X1dis, X1dis, r, c, mAdj, mImage, fDistanceFunc);
                    mImage(r,c) = minStepSize + mImage(r,c);
                    vX0 = fDistanceFunc.updateX0Add(vX0, vX1, minStepSize);
                end % end of inner for

            end % end of outer for       


        end % end of function
        

        function [minStepSize] = computeOptimalStepSize(obj, X0dis, X0X1dis, X1dis, basisRow, basisCol, mAdj, mImage, fObjective)
                    %
                    % Computes the distance as a function of stepSize, given the other
                    % traces are computed already.
                    %

                    a3 = fObjective.m_fImageReg.m_mA(basisRow, basisCol);
                    a2 = fObjective.m_fImageReg.m_mB(basisRow, basisCol);
                    a1 = fObjective.m_fImageReg.m_mC(basisRow, basisCol);
                    currM = mImage(basisRow, basisCol);
                    
                    quadB = 2 * X1dis + fObjective.m_regPara * (3 * a3 .* currM + a2);
                    quadC = -2 * X0X1dis - obj.m_l1ImageReg + fObjective.m_regPara * (3 * a3 .* currM^2 + 2 * a2 .* currM + a1);
                    
                    if quadB^2 - 12 * quadC < 0
                        % solve using line search
                        ftObj = @(t) obj.computeObjective(X0dis, X0X1dis, X1dis, basisRow, basisCol, mAdj, mImage, fObjective, t); 
                        [minStepSize, ~] = fminbnd(ftObj, -currM, 1 - currM);
                    else
                        sqPart = sqrt(quadB^2 - 12 * quadC);

                        stepSize1 = (-quadB - sqPart) / 6;
                        stepSize2 = (-quadB + sqPart) / 6;

                        % stepSize cannot be smaller than 0 or larger than 1
                        if stepSize1 < 0
                            stepSize1 = max(-currM, stepSize1);
                        else
                            stepSize1 = min(1-currM, stepSize1);
                        end
                        if stepSize2 < 0
                            stepSize2 = max(-currM, stepSize2);
                        else
                            stepSize2 = min(1-currM, stepSize2);
                        end                    

                        % now compare the two based on the objective
                        dis1 = obj.computeObjective(X0dis, X0X1dis, X1dis, basisRow, basisCol, mAdj, mImage, fObjective, stepSize1);
                        dis2 = obj.computeObjective(X0dis, X0X1dis, X1dis, basisRow, basisCol, mAdj, mImage, fObjective, stepSize2);

                        if dis1 < dis2
                            minStepSize = stepSize1;
                        else
                            minStepSize = stepSize2;
                        end
                    end
                    
        end % end of function        
        
        
        function [dis] = computeObjective(obj, X0dis, X0X1dis, X1dis, basisRow, basisCol, mAdj, mImage, fObjective, stepSize)
                    %
                    % Computes the distance as a function of stepSize, given the other
                    % traces are computed already.
                    %

                    dis = X0dis - 2 * stepSize * X0X1dis + stepSize * stepSize * X1dis - obj.m_l1ImageReg +...
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



