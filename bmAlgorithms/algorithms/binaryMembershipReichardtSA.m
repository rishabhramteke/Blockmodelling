function [mBestImage, mBestMembership, bestObjVal, totalIterNum] = binaryMembershipReichardtSA(mAdj, k, runNum, sDistanceFunc, varargin)
%
% Input:
% mAdj      - Adjacency matrix (n x n).
% k         - Number of positions to infer.
% convEpsilon   - Epsilon value to denote convergence.
% runNum    - Number of runs to perform.  Each run we start with a random
%             initial mImage.
% sUpdateFunc - The name of the update function to use.
%
% Output:
% mImage    - Image matrix (k x k).
% mMembership   - Membership matrix (n x k).
%
% Computes the hard approximation of A of blockmodel approximate mMem *
% mImage * mMem'.
% Uses Joerg Reichardt's null model formulation.
%
%
% @author: Jeffrey Chan, 2013
%


    % standard gamma
    gamma = 1.0;


    % initialisation for membership
    fMemInitFunc = InitRandomMem(false);
    

    % determine distance function
    switch sDistanceFunc
        case 'reichEqual'
            fDistanceFunc = ReichardtEqual(gamma);
            fHeatBathAlgor = @heatBathEqualNull;
        case 'reichDeg'
            fDistanceFunc = ReichardtDeg(gamma);
            fHeatBathAlgor = @heatBathDegNull;
        otherwise
            error('binaryMembershipReichardtSA:sDistanceFunc', 'Invalid distance function option');
    end

    totalIterNum = 0;


    % initial run
    if (runNum > 0)
        fprintf('run %d', 1);
        [mBestImage, mBestMembership, bestObjVal, bestMatApproxVal] = singleRun(mAdj, k, fDistanceFunc, fMemInitFunc, fHeatBathAlgor, varargin{:});
    else
        error('binaryMembershipBMAlgor:runNum', 'runNum is 0 or less, meaning algorithm was not executed for even one run');
    end


    % loop through remaining number of runs
    for r = 2 : runNum
        fprintf('run %d', r);

        [mCurrImage, mCurrMembership, currObjVal, bestMatApproxVal] = singleRun(mAdj, k, fDistanceFunc, fMemInitFunc, fHeatBathAlgor, varargin{:});

        % see if this is best solution and update as appropriate
        if (currObjVal < bestObjVal)
            bestObjVal = currObjVal;
            mBestImage = mCurrImage;
            mBestMembership = mCurrMembership;
        end
    end % end of outer for




end % end of function binarykMeansBMAlgor()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555

function [mImage, mMembership, currObjVal, currMatApproxVal] = singleRun(mAdj, k, fDistanceFunc, fMemInitFunc, fHeatBathAlgor, varargin)
%
% A single run of the algorithm.
%


    % parse arguments
    inParser = inputParser;
    inParser.KeepUnmatched = true;

    % add parameters and set default valuse
    % initial membership and image matrices
    addParameter(inParser, 'initialMem', NaN);
    
    parse(inParser, varargin{:});
        
    mInitMem = inParser.Results.initialMem;

    bInitMem = false;
    if ~isnan(mInitMem)
        bInitMem = true;
    end



    % initialise membership
    if ~bInitMem
        mMembership = fMemInitFunc.initMembership(mAdj, k);
    else
        mMembership = mInitMem;
    end    
    
    
    % options
%     gammaStart = 0.5;
%     gammaEnd = 2.0;
%     steps = 4;
%     repetitionNum = 10;
    
    stopTemperature = 1;
    coolFactor = 0.95;
    
    %initTemperature = (size(mAdj,1) + k * k) * 0.1;
    initTemperature = 1.0;
    
    % find optimal initial temperature
    disp('finding starting temperature');
    temperature = findStartTemperature(mAdj, mMembership, fDistanceFunc.m_gamma, initTemperature, fHeatBathAlgor);
    
    changeNum = 1;
    maxSweepNum = 50;
    while changeNum > 0 && (temperature / stopTemperature) > 1.0
        temperature = temperature * coolFactor;
        fprintf('current temperature = %4f', temperature);
        
        [acceptance, mMembership] = fHeatBathAlgor(mAdj, mMembership, temperature, fDistanceFunc.m_gamma, maxSweepNum);
        
        % check for stopping criteria
        if acceptance < (1 - 1/k)*0.01
            changeNum = 0;
        else
            changeNum = 1;
        end        
        
    end % end of while
    
       
       
    % obtain mImage
    mImage = computeImage(mAdj, mMembership, fDistanceFunc.m_gamma);
    currObjVal = fDistanceFunc.distance(mAdj, mMembership);
    
end % end of function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [mMembership] = gammaSweep(mAdj, mMembership, gammaStart, gammaStop, steps, repetitionNum, initTemperature, fHeatBathAlgor)
    %
    % Determine the optimal gamma?
    %

    stepSize = (gammaStop - gammaStart) / steps;
    posNum = size(mMembership, 2);
    
    for n = 1 : steps
        gamma = gammaStart + stepSize * (n-1);
        % initial temperature
        temperature = 0.5;
        
        acceptance = 0.5;
        maxSweepNum = 25;
        % from what i can infer, this is used to find the initial temperature
        while acceptance < (1 - 1/posNum) * 0.95
           temperature = temperature * 1.1;
           fHeatBathAlgor(mAdj, mMembership, temperature, gamma, maxSweepNum); 
        end
        
        % change maxswewep num for the rest of the iterations
        maxSweepNum = 50;
        for i = 1 : repetitionNum
           changeNum = 1;
           temperature = initTemperature;
           while changeNum > 0 && temperature > 0.01
              temperature = temperature * 0.99;
              [acceptance, mMembership] = fHeatBathAlgor(mAdj, mMembership, temperature, gamma, maxSweepNum);
              
              % check for stopping criteria
              if acceptance > (1 - 1/posNum)*0.01
                  changeNum = 1;
              else
                  changeNum = 0;
              end
           end
        end
        
        
        
    end % end of for

end



function [temperature] = findStartTemperature(mAdj, mMembership, gamma, initTemperature, fHeatBathAlgor)
    %
    % Finds the optimal starting temperature.
    %
    
    posNum = size(mMembership, 2);
    
    temperature = initTemperature;
    
    acceptance = 0.0;
    maxSweepNum = 50;

    
    % from what i can infer, this is used to find the initial temperature
    while acceptance < (1 - 1/posNum) * 0.95
       temperature = temperature * 1.1;
       mMembershipCopy = mMembership;
       acceptance = fHeatBathAlgor(mAdj, mMembershipCopy, temperature, gamma, maxSweepNum); 
    end    
    
    temperature = temperature * 1.1;
    
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [acceptance, mMembership] = heatBathEqualNull(mAdj, mMembershipIn, temperature, gamma, maxSweepNum)
%
% Generates a new neighbouring solution.
%
% cOptimValues - struct/cell structure of the current parameters of the Matlab
% SA algorithm.
%



%     for t = 1 : floor(cOptimValues.temperature) + 1
    sweepNum = 0;
    changesNum = 0;
    
    mMembership = mMembershipIn;
    
    memberNum = size(mMembership,1);
    
    while (sweepNum < maxSweepNum)
        sweepNum = sweepNum + 1;
       
        % compute the expected mImage matrix
        [mBlockSum, vElemNum, p] = computeExpectedImageMat(mAdj, mMembership);
        
        
        % perform a sweep
        % loop over all vertices in the sweep
        for v = 1 : memberNum
            % pick a random vertex
            pickedV = randi(memberNum, 1);
            
            vI = find(mMembership(pickedV,:) > 0);
            assert(size(vI,2) == 1 && size(vI,1) == 1);
            
            % current position
            oldPos = vI(1);
%             if size(oldPos,2) > 1 || size(oldPos,1) > 1
%                 pickedV
%                 mMembership
%             end
                
            beta = 1.0 / temperature;
                
            
            posNum = size(mMembership,2);
            mDeltaH = zeros(posNum, 1);
            % loop over the positions
            for newPos = 1 : posNum
                for q = 1 : posNum
                    % get the vertices of the positions
                    vPosVerts = find(mMembership(:, q) > 0);
                    
                    if size(abs(mBlockSum(newPos, q) + sum(mAdj(pickedV, vPosVerts)) - gamma * p * (vElemNum(newPos) + 1) * vElemNum(q)), 1) ~= 1
                        disp('4');
                        abs(mBlockSum(newPos, q) + sum(mAdj(pickedV, vPosVerts)) - gamma * p * (vElemNum(newPos) + 1) * vElemNum(q))
                    end
                    
                    if size(abs(mBlockSum(oldPos, q) + sum(mAdj(pickedV, vPosVerts)) - gamma * p * (vElemNum(oldPos) - 1) * vElemNum(q)), 1) ~= 1
                        disp('3');
                        abs(mBlockSum(oldPos, q) + sum(mAdj(pickedV, vPosVerts)) - gamma * p * (vElemNum(oldPos) - 1) * vElemNum(q))
                        oldPos
                        mBlockSum
                        mBlockSum(oldPos,q)
                        mAdj(pickedV, vPosVerts)
                        gamma
                        p
                        vElemNum(oldPos)
                        vElemNum(q)
                    end
                    
                    if size(abs(mBlockSum(newPos, q) - gamma *  p * vElemNum(newPos) * vElemNum(q)), 1) ~= 1
                        disp('2');
                        abs(mBlockSum(newPos, q) - gamma *  p * vElemNum(newPos) * vElemNum(q))
                    end
                    
                    if size(abs(mBlockSum(oldPos, q) - gamma * p * vElemNum(oldPos) * vElemNum(q))) ~= 1
                        disp('1');
                        abs(mBlockSum(oldPos, q) - gamma * p * vElemNum(oldPos) * vElemNum(q))
                    end
                    
                    
                    % out
                    mDeltaH(newPos) = mDeltaH(newPos) + abs(mBlockSum(oldPos, q) - gamma * p * vElemNum(oldPos) * vElemNum(q)) +...
                        abs(mBlockSum(newPos, q) - gamma *  p * vElemNum(newPos) * vElemNum(q)) -...
                        abs(mBlockSum(oldPos, q) + sum(mAdj(pickedV, vPosVerts)) - gamma * p * (vElemNum(oldPos) - 1) * vElemNum(q)) -...
                        abs(mBlockSum(newPos, q) + sum(mAdj(pickedV, vPosVerts)) - gamma * p * (vElemNum(newPos) + 1) * vElemNum(q));
                    % in
                    mDeltaH(newPos) = mDeltaH(newPos) + abs(mBlockSum(q, oldPos) - gamma * p * vElemNum(q) * vElemNum(oldPos)) +...
                        abs(mBlockSum(q, newPos) - gamma * p * vElemNum(q) * vElemNum(oldPos)) -...
                        abs(mBlockSum(q, oldPos) - sum(mAdj(vPosVerts, pickedV)) - gamma * p * vElemNum(q) * (vElemNum(oldPos) - 1) ) -...
                        abs(mBlockSum(q, newPos) - sum(mAdj(vPosVerts, pickedV)) - gamma * p * vElemNum(q) * (vElemNum(newPos) + 1) );
                    
                    % normalise
                    mDeltaH(newPos) = mDeltaH(newPos) / (size(mAdj,1))^2;
                end
                
            end
            
            % pick the appropriate position 
            % compute sum of exp^H
            totalH = 0;
            for p = 1 : posNum
                totalH = totalH + exp(-beta * mDeltaH(p,1));
            end
            
            vProbH = zeros(posNum, 1);
            for p = 1 : posNum
               vProbH(p) = exp(-beta * mDeltaH(p,1)) / totalH; 
            end
            
            % generate unif random number in [0,1]
            r = rand(1);
            
            % see which thing prob range r falls into
            currSum = 0.0;
            choosenPos = -1;
            for p = 1 : posNum
                currSum = currSum + vProbH(p);
                if r <= currSum
                    choosenPos = p;
                    break; 
                end
            end
            
            if (choosenPos <= 0)
                disp(sprintf('r = %d', r));
                disp(sprintf('currSum = %4f', currSum));
                vProbH
                mDeltaH
                beta
                totalH
            end
            assert(choosenPos > 0);
                 
            % update the positions
            if choosenPos ~= oldPos
               changesNum = changesNum + 1;
               
               % must update blocksum before mMembership is updated
               mBlockSum = updateBlockSum(mBlockSum, mAdj, mMembership, pickedV, oldPos, choosenPos);
               
               if size(find(mMembership(pickedV, :) > 0), 2) > 1 || size(find(mMembership(pickedV, :) > 0), 1) > 1
                   disp('before');
                   pickedV
                   oldPos
                   choosenPos
                   mMembership(pickedV,:)
               end               
               
               
               mMembership(pickedV, oldPos) = 0;
               mMembership(pickedV, choosenPos) = 1;
               
               if size(find(mMembership(pickedV, :) > 0), 2) > 1 || size(find(mMembership(pickedV, :) > 0), 1) > 1
                   disp('after');
                   pickedV
                   oldPos
                   choosenPos
                   mMembership(pickedV,:)
               end
               assert(size(find(mMembership(pickedV, :) > 0), 2) == 1 && size(find(mMembership(pickedV, :) > 0), 1) == 1);
               
               vElemNum(oldPos) = vElemNum(oldPos) - 1;
               vElemNum(choosenPos) = vElemNum(choosenPos) + 1;

               
               

            end
        end % end of for over number of vertices
       
%     end % end of for
    end % end of while

    
    
    acceptance = changesNum / (memberNum * sweepNum);
end % end of function



function [acceptance, mMembership] = heatBathDegNull(mAdj, mMembershipIn, temperature, gamma, maxSweepNum)
%
% Generates a new neighbouring solution.
%
% cOptimValues - struct/cell structure of the current parameters of the Matlab
% SA algorithm.
%



%     for t = 1 : floor(cOptimValues.temperature) + 1
    sweepNum = 0;
    changesNum = 0;
    
    mMembership = mMembershipIn;
    
    memberNum = size(mMembership,1);
    
    while (sweepNum < maxSweepNum)
        sweepNum = sweepNum + 1;
       
        % compute the expected mImage matrix
        [mBlockSum] = computeExpectedImageMat(mAdj, mMembership);
        
        % compute the total in and out degree of each position
        [vOutDegSum, vInDegSum] = computeTotalDeg(mAdj, mMembership);
        
        
        % perform a sweep
        % loop over all vertices in the sweep
        for v = 1 : memberNum
            % pick a random vertex
            pickedV = randi(memberNum, 1);
            
            vI = find(mMembership(pickedV,:) > 0);
            assert(size(vI,2) == 1 && size(vI,1) == 1);
            
            % current position
            oldPos = vI(1);
%             if size(oldPos,2) > 1 || size(oldPos,1) > 1
%                 pickedV
%                 mMembership
%             end
                
            beta = 1.0 / temperature;
                
            
            posNum = size(mMembership,2);
            mDeltaH = zeros(posNum, 1);
            % loop over the positions
            for newPos = 1 : posNum
                for q = 1 : posNum
                    % get the vertices of the positions
                    vPosVerts = find(mMembership(:, q) > 0);
                                        
                    outVtoQ = sum(mAdj(pickedV, vPosVerts));
                    inQtoV = sum(mAdj(vPosVerts, pickedV));
                    
                    % out
                    mDeltaH(newPos) = mDeltaH(newPos) + abs(mBlockSum(oldPos, q) - gamma * vOutDegSum(oldPos) * vInDegSum(q)) +...
                        abs(mBlockSum(newPos, q) - gamma *  vOutDegSum(newPos) * vInDegSum(q)) -...
                        abs(mBlockSum(oldPos, q) - outVtoQ - gamma * (vOutDegSum(oldPos) - outVtoQ) * vInDegSum(q)) -...
                        abs(mBlockSum(newPos, q) + outVtoQ - gamma * (vOutDegSum(newPos) + outVtoQ) * vInDegSum(q));
                    % in
                    mDeltaH(newPos) = mDeltaH(newPos) + abs(mBlockSum(q, oldPos) - gamma * vOutDegSum(q) * vInDegSum(oldPos)) +...
                        abs(mBlockSum(q, newPos) - gamma * vOutDegSum(q) * vInDegSum(oldPos)) -...
                        abs(mBlockSum(q, oldPos) - inQtoV - gamma * vOutDegSum(q) * (vInDegSum(oldPos) - inQtoV) ) -...
                        abs(mBlockSum(q, newPos) + inQtoV - gamma * vOutDegSum(q) * (vInDegSum(newPos) + inQtoV) );
                    
                    % normalise
                    mDeltaH(newPos) = mDeltaH(newPos) / (size(mAdj,1))^2;
                end
                
            end
            
            % pick the appropriate position 
            % compute sum of exp^H
            totalH = 0;
            for p = 1 : posNum
                totalH = totalH + exp(-beta * mDeltaH(p,1));
            end
            
            vProbH = zeros(posNum, 1);
            for p = 1 : posNum
               vProbH(p) = exp(-beta * mDeltaH(p,1)) / totalH; 
            end
            
            % generate unif random number in [0,1]
            r = rand(1);
            
            % see which thing prob range r falls into
            currSum = 0.0;
            choosenPos = -1;
            for p = 1 : posNum
                currSum = currSum + vProbH(p);
                if r <= currSum
                    choosenPos = p;
                    break; 
                end
            end
            
            if (choosenPos <= 0)
                disp(sprintf('r = %d', r));
                disp(sprintf('currSum = %4f', currSum));
                vProbH
                mDeltaH
                beta
                totalH
            end
            assert(choosenPos > 0);
                 
            % update the positions
            if choosenPos ~= oldPos
               changesNum = changesNum + 1;
               
               % must update blocksum before mMembership is updated
               mBlockSum = updateBlockSum(mBlockSum, mAdj, mMembership, pickedV, oldPos, choosenPos);
               
               % update degree
               [vOutDegSum, vInDegSum] = updateTotalDeg(vOutDegSum, vInDegSum, mAdj, pickedV, oldPos, choosenPos);
               
               
               if size(find(mMembership(pickedV, :) > 0), 2) > 1 || size(find(mMembership(pickedV, :) > 0), 1) > 1
                   disp('before');
                   pickedV
                   oldPos
                   choosenPos
                   mMembership(pickedV,:)
               end               
               
               
               mMembership(pickedV, oldPos) = 0;
               mMembership(pickedV, choosenPos) = 1;
               
               if size(find(mMembership(pickedV, :) > 0), 2) > 1 || size(find(mMembership(pickedV, :) > 0), 1) > 1
                   disp('after');
                   pickedV
                   oldPos
                   choosenPos
                   mMembership(pickedV,:)
               end
               assert(size(find(mMembership(pickedV, :) > 0), 2) == 1 && size(find(mMembership(pickedV, :) > 0), 1) == 1);
                              

            end
        end % end of for over number of vertices
       
%     end % end of for
    end % end of while

    
    
    acceptance = changesNum / (memberNum * sweepNum);
end % end of function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55


function [mBlockSum, vElemNum, edgeProb] = computeExpectedImageMat(mAdj, mMembership)
%
% Computes the expected image matrix according to the equal edge probability
% null model.

    posNum = size(mMembership, 2);
    % update m_p object
    edgeProb = sum(sum(mAdj)) / (size(mAdj,1))^2;
        
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
%     mExpected = zeros(posNum, posNum);
%     for r = 1 : posNum
%         for c = 1 : posNum
%             mExpected(r,c) = edgeProb * vElemNum(r) * vElemNum(c);
%         end
%     end
    
end % end of function



function [vOutDegSum, vInDegSum] = computeTotalDeg(mAdj, mMembership)
%
% Computes the expected image matrix according to the equal edge probability
% null model.

    posNum = size(mMembership, 2);
    
    vOutDegSum = zeros(1, posNum);
    vInDegSum = zeros(1, posNum);
    
    % loop through each edge
    [vR, vC] = find(mAdj);
    
    % check if vR is a row or column vector
    [nr, nc] = size(vR);
    nonZeroNum = nr;
    if nr < nc
        nonZeroNum = nc;
    end
    
    for e = 1 : nonZeroNum
        % find the position of the incident vertices
        srcPos = find(mMembership(vR(e), :) > 0);
        assert(size(srcPos,1) == 1);
        desPos = find(mMembership(vC(e), :) > 0);
        assert(size(desPos,1) == 1);
        
        % update the deg sum vectors
        vOutDegSum(srcPos) = vOutDegSum(srcPos) + mAdj(vR(e), vC(e));
        vInDegSum(desPos) = vInDegSum(srcPos) + mAdj(vR(e), vC(e));
    end

    
end % end of function




function [mBlockSum] = updateBlockSum(mBlockSum, mAdj, mMembership, v, oldPos, newPos)
%
% Update the block sum for vertex pickedV moving from position 'oldPos' to
% 'newPos'
%
    assert(oldPos ~= newPos);

    [vRow] = find(mAdj(v, :) > 0);
    [vCol] = find(mAdj(:, v) > 0);
    
    % row (TODO check vRow and vCol are transversed in the right dimensions)
    for oe = 1 : size(vRow, 2)
        % find the position that vRow(oe) is currently in
        vPos = find(mMembership(vRow(oe),:) > 0);
        assert(size(vPos,1) == 1);
        
        mBlockSum(oldPos, vPos(1)) = mBlockSum(oldPos, vPos(1)) - mAdj(v, vRow(oe));
        mBlockSum(newPos, vPos(1)) = mBlockSum(newPos, vPos(1)) + mAdj(v, vRow(oe));
    end
        
    % column
    for ie = 1 : size(vCol, 1)
        % find the position that vRow(ie) is currently in
        vPos = find(mMembership(vCol(ie),:) > 0);
        assert(size(vPos,1) == 1);
        
        mBlockSum(vPos(1), oldPos) = mBlockSum(vPos(1), oldPos) - mAdj(vCol(ie), v);
        mBlockSum(vPos(1), newPos) = mBlockSum(vPos(1), newPos) + mAdj(vCol(ie), v);
    end
    
    % see if we subtracted twice
    if mAdj(v,v) > 0
        mBlockSum(oldPos, oldPos) = mBlockSum(oldPos, oldPos) + mAdj(v,v);
        mBlockSum(newPos, newPos) = mBlockSum(newPos, newPos) - mAdj(v,v);
    end
    


end %end of function




function [vOutDegSum, vInDegSum] = updateTotalDeg(vOutDegSum, vInDegSum, mAdj, pickedV, oldPos, newPos)
%
% Update the total degree of positions based on pickedV moving from oldPos to
% newPos.
%
    outDegSum = sum(mAdj(pickedV, :));
    inDegSum = sum(mAdj(:, pickedV));
    
    vOutDegSum(oldPos) = vOutDegSum(oldPos) - outDegSum;
    vInDegSum(oldPos) = vInDegSum(oldPos) - inDegSum;
    
    vOutDegSum(newPos) = vOutDegSum(newPos) - outDegSum;
    vInDegSum(newPos) = vInDegSum(newPos) - inDegSum;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [dis] = objectiveFunc(vSolution, mAdj, fObjFunc, membershipNum, posNum)
%
% Objective function for the SA to optimise.
%
    % cSolution{2} = mImage, cSolution{1} = mMembership
    dis = fObjFunc.quality(mAdj, reshape(vSolution(1:membershipNum*posNum), membershipNum, posNum));
    
    
end % end of function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vSolution] = linearise(mMembership)
%
% Linearise the matrix into one vector.  Columnwise linearisation.
%
    vSolution = mMembership(:);
    
end % end of function



function [mImage] = computeImage(mAdj, mMembership, gamma)
% 
% Computes the image matrix.
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
            
    mImage = zeros(posNum, posNum);
    for r = 1 : posNum
        for s = 1 : posNum
            if mBlockSum(r,s) - gamma * mExpected(r,s) > 0
                mImage(r,s) = 1;
            else
                mImage(r,s) = 0;
            end
        end
    end

end % end of function


function val = getValue(vSolution, bMembership, row, col, membershipNum, posNum)
%
% Returns the element at matrix(row,col), where matrix is mMembership if
% bMembership = true or mImage otherwise.
%
    if bMembership
       val = vSolution(row + (col - 1) * membershipNum); 
    else
       val = vSolution(membershipNum * posNum + row + (col-1) * posNum); 
    end
end % end of function


function [vSolution] = setValue(vSolution, bMembership, row, col, val, membershipNum, posNum)
%
% Sets the element at matrix(row,col), where matrix is mMembership if
% bMembership = true or mImage otherwise.
%
    if bMembership
       vSolution(row + (col - 1) * membershipNum) = val;
    else
       vSolution(membershipNum * posNum + row + (col-1) * posNum) = val;
    end
end % end of function


function [vRow] = getMembershipRow(vSolution, row, membershipNum, posNum)
%
% Returns the row 'row' of mMembership.
%
    
    vRow = vSolution([row : membershipNum : row + (posNum - 1) * membershipNum]); 

end
