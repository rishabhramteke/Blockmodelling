function [mBestImage, mBestMembership, bestObjVal, bestMatApproxVal, totalIterNum] = binaryMembershipSA(mAdj, k, runNum, sDistanceFunc)
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
% Uses simulated annealing.
%
%
% @author: Jeffrey Chan, 2013
%


totalIterNum = 0;

fInitFunc = @initRandomImage;

% determine distance function
switch sDistanceFunc
    case 'euclidean'
        fDistanceFunc = EucTriFact;
        fMatApproxDisFunc = EucTriFact;
    case 'bmEuclidean'
        fDistanceFunc = EucTriFactBM;
        fMatApproxDisFunc = EucTriFact;
    otherwise
        error('binaryMembershipBMAlgor:sDistanceFunc', 'Invalid distance function option');
end


% initial run
if (runNum > 0)
    disp(sprintf('run %d', 1));
    [mBestImage, mBestMembership, bestObjVal, bestMatApproxVal] = singleRun(mAdj, k, fInitFunc, fDistanceFunc, fMatApproxDisFunc);
else
    error('binaryMembershipBMAlgor:runNum', 'runNum is 0 or less, meaning algorithm was not executed for even one run');
end


% loop through remaining number of runs
for r = 2 : runNum
    disp(sprintf('run %d', r));
    
    [mCurrImage, mCurrMembership, currObjVal, bestMatApproxVal] = singleRun(mAdj, k, fInitFunc, fDistanceFunc, fMatApproxDisFunc);
    
    % see if this is best solution and update as appropriate
    if (currObjVal < bestObjVal)
        bestObjVal = currObjVal;
        mBestImage = mCurrImage;
        mBestMembership = mCurrMembership;
    end
end % end of outer for




end % end of function binarykMeansBMAlgor()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555

function [mImage, mMembership, currObjVal, currMatApproxVal] = singleRun(mAdj, k, fInitFunc, fDistanceFunc, fMatApproxDisFunc)
%
% A single run of the algorithm.
%

% MODIFICATION: We should update mImage first, given random initial
% membership, as we might get convergence problems where the mImage forces
% one or more positions to lose all its members


    % initial random membership matrix
    mMembership = randomInitMembership(size(mAdj,1), k);

    % initialise mImage
    mImage = fInitFunc(mAdj, mMembership);
    
    % options
    fNeighbour = @(optim, data) genNewSolution(optim, data, size(mAdj,1), k);
    saOptions = saoptimset('DataType', 'custom', 'AnnealingFcn', fNeighbour,...
        'InitialTemperature', (size(mAdj,1) + k * k) * 0.1,...
        'TemperatureFcn', @temperatureexp,...
        'ReannealInterval', 50);
    % initial temperature is the number of vetices + number of entries in the
    % image matrix and 10%
    % temperature is changed after 50 iterations
    % use default temperature cooling scheme
    
    
    % intialise (random) solution (linearise it)
    vInitSolution = linearise(mMembership, mImage);
    %
    fFitness = @(vSolution) objectiveFunc(vSolution, mAdj, fDistanceFunc, size(mAdj,1), k);
    
    
    [vNewSolution, currObjVal] = simulannealbnd(fFitness, vInitSolution, [], [], saOptions);
    
    mImage = reshape(vNewSolution(size(mAdj,1)*k+1:end), k, k);
    mMembership = reshape(vNewSolution(1:size(mAdj,1)*k), size(mAdj,1), k);
    % harden mImage
%     mImage = full(mImage);
    currMatApproxVal = fMatApproxDisFunc.distance(mAdj, mImage, mMembership);
    


end % end of function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [vSolution] = genNewSolution(cOptimValues, problemData, memberNum, posNum)
%
% Generates a new neighbouring solution.
%
% cOptimValues - struct/cell structure of the current parameters of the Matlab
% SA algorithm.
%

    % parameters
    imagePropPara = 0.5;

    for t = 1 : floor(cOptimValues.temperature) + 1
        % choose either to optimise membership or image
        if (rand(1) < imagePropPara)
            % permute image
            vSolution = changeImage(cOptimValues.x, memberNum, posNum);
        else
            % permute membership
            vSolution = changeMembership(cOptimValues.x, memberNum, posNum);
        end
    end % end of for
    

end % end of function



function [vSolution] = changeImage(vSolution, memberNum, posNum)
%
% Randomly update image matrix.
% Hard image.
%
    % generate the entry to change
    row = randi(posNum, 1);
    col = randi(posNum, 1);
    
    % discrete
    if getValue(vSolution, false, row, col, memberNum, posNum) > 0
        % set to 0
        vSolution = setValue(vSolution, false, row, col, 0, memberNum, posNum);
    else
        % set to 1
        vSolution = setValue(vSolution, false, row, col, 1, memberNum, posNum);
    end

end % end of function


function [vSolution] = changeMembership(vSolution, memberNum, posNum)
%
% Randomly update membership matrix.
%

    % generate the entry to change
    row = randi(memberNum, 1);
    newPos = randi(posNum, 1);
    
    % discrete
    currPos = find(getMembershipRow(vSolution, row, memberNum, posNum) > 0);

    vSolution = setValue(vSolution, true, row, currPos, 0, memberNum, posNum);
    vSolution = setValue(vSolution, true, row, newPos, 1, memberNum, posNum);
  
end % end of function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [dis] = objectiveFunc(vSolution, mAdj, fObjFunc, membershipNum, posNum)
%
% Objective function for the SA to optimise.
%
    % cSolution{2} = mImage, cSolution{1} = mMembership
    dis = fObjFunc.distance(mAdj, reshape(vSolution(membershipNum*posNum+1:end), posNum, posNum), reshape(vSolution(1:membershipNum*posNum), membershipNum, posNum));
    
end % end of function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vSolution] = linearise(mMembership, mImage)
%
% Linearise the two matrices into one vector.  Columnwise linearisation.
%
    vSolution = cat(1, mMembership(:), mImage(:));
    
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
