classdef ReichardtDeg
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
        
        function obj = ReichardtDeg(gamma)
            obj.m_gamma = gamma;            
        end % end of constructor()
        
        
        function [qual] = distance(obj, mAdj, mMembership)
        %
        % Computes the quality function.
        %
        
            posNum = size(mMembership, 2);
            vOutDegSum = zeros(1, posNum);
            vInDegSum = zeros(1, posNum);

            
            % find non-zero entries for adjacency matrix
            [vNZRows, vNZCols] = find(mAdj);
            
            % compute m_{rs} (sum of edge weights for block rs
            mBlockSum = zeros(posNum, posNum);
            % go through each non zero entry
            for e = 1 : size(vNZRows,1)
                rowPos = find(mMembership(vNZRows(e),:));
                colPos = find(mMembership(vNZCols(e),:));
                
                mBlockSum(rowPos, colPos) = mBlockSum(rowPos, colPos) + mAdj(vNZRows(e), vNZCols(e));
                
                % update the deg sum vectors
                vOutDegSum(rowPos) = vOutDegSum(rowPos) + mAdj(vNZRows(e), vNZCols(e));
                vInDegSum(colPos) = vInDegSum(colPos) + mAdj(vNZRows(e), vNZCols(e));                
            end
            
            % compute the sum of edges
            edgeSum = 1 / sum(sum(mAdj));
            
            % compute the expected matrix
            mExpected = zeros(posNum, posNum);
            for r = 1 : posNum
                for c = 1 : posNum
                    mExpected(r,c) = edgeSum * vOutDegSum(r) * vInDegSum(c);
                end
            end
                     
            qual = 0.5 * sum(sum(abs(mBlockSum - obj.m_gamma * mExpected)));
        
        end % end of distance()
        
        

    end 

end % end of class
