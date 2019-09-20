classdef EucImageDistance
    %
    % Compate between input image and reference image.
    %
    % @author: Jeffrey Chan, 2014
    %
    
    properties
        m_mRefImage;
        m_mRefPos;
    end
    
    
    methods
        function [obj] = EucImageDistance(mRefImage, mRefPos)
            obj.m_mRefImage = mRefImage;
            obj.m_mRefPos = mRefPos;
        end
        
        
        function [dis] = distance(obj, mImage, mPos)
            % align the positions
            % for each permutation of the positions, we want to find the closest
            % permutation.  This permutation is then used to rearrange the image
            % matrix and then compared with ideal.
            
            mPerm = perms(1:size(obj.m_mRefPos,2));
            mPerm = mPerm(:,end:-1:1);
            minPermDis = sum(sum(abs(mPos(:,mPerm(1,:)) - obj.m_mRefPos)));
            minPermIndex = 1;
            for i = 2 : size(mPerm,1)
                permDis = sum(sum(abs(mPos(:,mPerm(i,:)) - obj.m_mRefPos)));
                if permDis < minPermDis
                    minPermDis = permDis;
                    minPermIndex = i;
                end
            end
            
            % realign mImage first
            mDiff = mImage(mPerm(minPermIndex,:), mPerm(minPermIndex,:)) - obj.m_mRefImage;
            dis = trace(mDiff*mDiff');
        end
    end
    
    
end % end of class