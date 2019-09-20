classdef ReichardtEqual
    %
    % Computes the Reichardt objective with a null model that assumes equal
    % probability for all edges.
    %
    %
    % @author: Jeffrey Chan, 2013
    %    

    properties
       
       % gamma parameter to weight how much mismatches are penalised in
       % comparison to matches
       m_gamma = 1
        
    end
    
    
    methods
        
        function obj = ReichardtEqual(gamma)
            obj.m_gamma = gamma;            
        end % end of constructor()
        
        
        function [qual] = distance(obj, mAdj, mMembership)
        %
        % Computes the quality function.
        %
        
            posNum = size(mMembership, 2);
            % update m_p object
            p = sum(sum(mAdj)) / (size(mAdj,1))^2;
        
            % sum up the number of elements in position
            vElemNum = zeros(posNum,1);
            for p = 1 : posNum
               vElemNum(p) = size(find(mMembership(:,p) > 0),1);
            end
            
            % find non-zero entries for adjacency matrix
            [vNZRows, vNZCols] = find(mAdj);
            
            % compute m_{rs} (sum of edge weights for block rs
            mBlockSum = zeros(posNum, posNum);
            % go through each non zero entry
            for e = 1 : size(vNZRows,1)
                rowPos = find(mMembership(vNZRows(e),:));
                colPos = find(mMembership(vNZCols(e),:));
                
                mBlockSum(rowPos, colPos) = mBlockSum(rowPos, colPos) + mAdj(vNZRows(e), vNZCols(e));
            end
            
            % compute the expected matrix
            mExpected = zeros(posNum, posNum);
            for r = 1 : posNum
                for c = 1 : posNum
                    mExpected(r,c) = p * vElemNum(r) * vElemNum(c);
                end
            end
            
            qual = 0.5 * sum(sum(abs(mBlockSum - obj.m_gamma * mExpected)));
        
        end % end of distance()
        
        

    end 

end % end of class
