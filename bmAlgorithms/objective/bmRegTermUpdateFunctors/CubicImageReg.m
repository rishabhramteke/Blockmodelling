classdef CubicImageReg < handle
    %
    % Cubic image regularisation function object.
    %
    %
    % @author: Jeffrey Chan, 2014
    %    
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        % maximum penalty value
        m_mMu = 0.5;
        % cubic parameters (ax^3 + bx^2 + cx + d)
        m_mA;
        m_mB;
        m_mC;
        m_mD;
    end % end of properties
    
    
    methods
            
        function obj = CubicImageReg(varargin)
            if length(varargin) == 1
                obj.m_mMu = varargin{1};
            end             
            
            % initialise structures
            obj.m_mA = zeros(size(obj.m_mMu,1), size(obj.m_mMu,2));
            obj.m_mB = zeros(size(obj.m_mMu,1), size(obj.m_mMu,2));
            obj.m_mC = zeros(size(obj.m_mMu,1), size(obj.m_mMu,2));
            obj.m_mD = zeros(size(obj.m_mMu,1), size(obj.m_mMu,2));
            
            % compute m_a, m_b, m_c, and m_d
            for r = 1 : size(obj.m_mMu,1)
                for c = 1 : size(obj.m_mMu,2)
                    mu = obj.m_mMu(r,c);
                    vX = [0, mu / 2, mu, mu + (1-mu)/2, 1];
                    vY = [0, obj.absoluteFunction(mu / 2, mu), 1, obj.absoluteFunction(mu + (1-obj.m_mMu(r,c)) / 2, mu), 0];

                    [vP] = polyfit(vX, vY, 3);
                    obj.m_mA(r,c) = vP(1);
                    obj.m_mB(r,c) = vP(2);
                    obj.m_mC(r,c) = vP(3);
                    obj.m_mD(r,c) = vP(4);
                end
            end
            
        end % end of function
        
        

            
        
        function [y] = absoluteFunction(~, m, mu)
            %
            % The absolute function we are trying to approximate
            %
            
            if m <= mu
                y = m / mu;
            else
                y = (1 - m) / (1 - mu);
            end
            
        end % end of function
        
                
        
        function dis = computeImageRegObj(obj, mImage)
            %
            % Compute the contribution of the regularisation term to the overall
            % objective.
            %
            
            mSq = mImage .* mImage;
            mCub = mImage .* mSq;
            dis = sum(sum(obj.m_mA .* mCub + obj.m_mB .* mSq + obj.m_mC .* mImage + obj.m_mD));
        end % end of function
        
        
        function [mImageRegGrad] = computeImageRegGrad(obj, mImage)
            %
            % Compute the gradient of the image regularisation term.
            % Does not factor in the regularisation weight.
            %
                        
            mImageRegGrad = 3 * obj.m_mA .* mImage .* mImage + 2 * obj.m_mB .* mImage + obj.m_mC;
        end % end of function()
    
        
        
        function dis = computeImageRegCoord(obj, mImage, stepSize, basisRow, basisCol)
            %
            % Compute the image regularisation objective for coordinate descent.
            %
            
            stepSq = stepSize * stepSize;
            stepCub = stepSize * stepSq;
                
            dis = stepCub + stepSq * (3 * obj.m_mA(basisRow, basisCol) + mImage(basisRow, basisCol) + obj.m_mB(basisRow, basisCol)) +...
                stepSize * (3 * obj.m_mA(basisRow, basisCol) * mImage(basisRow, basisCol)^2 + 2 * obj.m_mB(basisRow, basisCol) * mImage(basisRow, basisCol) + obj.m_mC(basisRow, basisCol)) +...
                (obj.m_mA(basisRow, basisCol) * mImage(basisRow, basisCol)^3 + obj.m_mB(basisRow, basisCol) * mImage(basisRow, basisCol)^2 + obj.m_mC(basisRow, basisCol) * mImage(basisRow, basisCol) + obj.m_mD(basisRow, basisCol));
            
        end % end of function
        
        
        
        
    end % end of methods
    
end % end of class