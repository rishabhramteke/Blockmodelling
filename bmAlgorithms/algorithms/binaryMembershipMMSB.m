function [mBestImage, mBestMembership, bestObjVal, totalIterNum] = binaryMembershipMMSB(mAdj, k, runNum, epsilonConv)
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
% Computes mixed membership stochastic blockmodel.
%


    fDistanceFunc = EucTriFactBM;


    totalIterNum = 0;

    % initial run
    if (runNum > 0)
        fprintf('run %d', 1);
        [mBestImage, mBestMembership, bestObjVal, iterNum] = singleRun(mAdj, k, epsilonConv, fDistanceFunc);
        totalIterNum = totalIterNum + iterNum;
    else
        error('binaryMembershipMMSB:runNum', 'runNum is 0 or less, meaning algorithm was not executed for even one run');
    end


    % loop through remaining number of runs
    for r = 2 : runNum
        fprintf('run %d', r);

        [mCurrImage, mCurrMembership, currObjVal, iterNum] = singleRun(mAdj, k, epsilonConv, fDistanceFunc);
        totalIterNum = totalIterNum + iterNum;

        % see if this is best solution and update as appropriate
        if (currObjVal < bestObjVal)
            bestObjVal = currObjVal;
            mBestImage = mCurrImage;
            mBestMembership = mCurrMembership;
        end
    end % end of outer for

end % end of function binarykMeansBMAlgor()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555


function [mImage, mMembership, currObjVal, iterNum] = singleRun(mAdj, k, epsilonConv, fDistanceFunc)
%
% A single run of the algorithm.
%

    vertNum = size(mAdj,1);
    iterNum = 1;
    
    
    % initialise parameters
    paraSparsity = sum(sum(mAdj)) / (vertNum * vertNum);
    mParaMulinomIn = ones(vertNum, vertNum, k) * (2 * vertNum / k);
    mParaMulinomOut = ones(vertNum, vertNum, k) * (2 * vertNum / k);
    mParaDirichlet = rand(vertNum, k);
    mParaDirichletAlpha = rand(vertNum, k);
    
    % M step
    mParaDirichletAlpha = updateDirichletAlpha(mParaDirichletAlpha, mParaDirichlet);
    mImage = updateImage(mAdj, mParaMulinomOut, mParaMulinomIn, k, paraSparsity);
    
    % main loop
    while epsilonConv
        for p = 1 : vertNum
            for q = 1 : vertNum
                % E step
                [mParaMultinomOut, mParaMultinomIn] =...
                    varationalInference(mAdj(p,q), vertNum, k, mParaDirichlet(p,:), mParaDirichlet(q,:), mParaMultinomIn, mParaMultinomIn, mImage, epsilonConv);
                mParaDirichlet = updateDirichletVar(mParaDirichletAlpha, mParaMultinomOut, mParaMultinomIn);
                % M step
                mParaDirichletAlpha = updateDirichletAlpha(mParaDirichletAlpha, mParaDirichlet);
                mImage = updateImage(mAdj, cmParaMulinomOut, cmParaMulinomIn, k, paraSparsity);                
            end
        end
        
        iterNum = iterNum + 1;
        cmParaMulinomOut{iterNum} = mParaMultinomOut;
        cmParaMulinomIn{iterNum} = mParaMultinomIn;
    end % end of while loop
    
    
    % get pi membership vectors
    
    mMembership = zeros(vertNum, k);
    for v = 1 : vertNum
        mMembership(v,:) = drchrnd(mParaDirichlet(v,:),1);
    end
    currObjVal = fDistanceFunc.distance(mAdj, mImage, mMembership);
    
end % end of function


function mImage = updateImage(mAdj, mParaMultinomOut, mParaMultinomIn, k, paraSparsity)
%
% M step for updating the image matrix variable
%

    
    mImage = zeros(k, k);

    for g = 1 : k
        for h = 1 : k
            mTemp = mParaMultinomOut(:,:,g) .* mParaMultinomIn(:,:,h);
            mImage(g,h) = mImage(g,h) + sum(sum(mAdj .* mTemp)) / ((1-paraSparsity) * sum(sum(mTemp)));   
        end
    end

end % end of function



function [vParaMultinomOut, vParaMultinomIn] =...
    varationalInference(edge, vertNum, k, vParaDirichletSrc, vParaDirichletTar, mParaMultinomOut, mParaMultinomIn, mImage, epsilonConv)

    vParaMultinomOut = ones(1,k);
    vParaMultinomIn = ones(1,k);

    vExpectedMembershipSrc = expectedMembership(mParaDirichletSrc);
    vExpectedMembershipTar = expectedMembership(vParaDirichletTar);
    
    while epsilonConv
        
        for g = 1 : k
            vParaMultinomOut(g) = exp(vExpectedMembershipSrc(g)) * prod(((mImage(g,:) .^ edge) .* ((1 - mImage(g,:) .^ (1 - edge)))) .^ mParaMultinomOut);
            vParaMultinomIn(g) = exp(vExpectedMembershipTar(g)) * prod(((mImage(:,g) .^ edge) .* ((1 - mImage(:,g) .^ (1 - edge)))) .^ mParaMultinomIn);
            % normalise to sum to 1
            vParaMultinomOut = vParaMultinomOut ./ sum(vParaMultinomOut);
            vParaMultinomIn = vParaMultinomIn ./ sum(vParaMultinomIn);
        end
        
        
    end
end % end of function



function vExpectedMembership = expectedMembership(vParaDirichlet)

    vExpectedMembership = psi(vParaDirichlet) - psi(sum(vParaDirichlet));
end % end of function


function mParaDirichlet = updateDirichletVar(mParaDirichletAlpha, mParaMultinomOut, mParaMultinomIn)
%
% Update Dirichlet variables.
%
    mParaDirichlet = repmat(mParaDirichletAlpha, size(mParaMultinomOut,1), 1) +...
        sum(mParaMultinomOut,2) + mParaMultinomIn(mParaMultinomOut,2);


end % end of function



function vParaDirichletAlpha = updateDirichletAlpha(mParaDirichletAlpha, mParaDirichlet)
%
% Update M step for the hyper parameters of the dirichlet distribution on
% membership variables.
%
    
    vertNum = size(mParaDirichlet, 1);
    
    options = optimoptions(@fminunc, 'GradObj', 'on', 'Hessian', 'on', 'Algorithm', 'quasi-newton', 'Display', 'off');
    for v = 1 : vertNum
        v
        f = @(x)dirichletParaFun(x, mParaDirichlet, v);
        vX0 = mParaDirichletAlpha(v,:);
        [vMinAlpha, fval, exitflag, output] = fminunc(f, vX0, options);
        mParaDirichletAlpha(v,:) = vMinAlpha;
    end
    
    
end % end of function


function [objVal, vGrad, mHessian] = dirichletParaFun(vParaDirichletAlpha, mParaDirichlet, v)
%
% Function of lower bound for Likelihood estimator.
%
    vertNum = size(mParaDirichlet,1);
    vParaDirichletAlpha
    
    lastTerm = sum( (vParaDirichletAlpha - 1) .* (psi(mParaDirichlet(v,:)) - repmat(psi(sum(mParaDirichlet(v,:))), 1, size(mParaDirichlet,2))) );
    objVal = log(gamma(sum(vParaDirichletAlpha))) - sum(log(gamma(vParaDirichletAlpha))) + lastTerm;
        
    vGrad = vertNum * (psi(sum(vParaDirichletAlpha)) - psi(vParaDirichletAlpha)) + sum(mParaDirichlet - repmat(psi(sum(mParaDirichlet,2)), 1, size(mParaDirichlet,2)), 1);
    
    mHessian = repmat(-vertNum * psi(1, sum(vParaDirichletAlpha)), length(vParaDirichletAlpha), length(vParaDirichletAlpha)) + vertNum * diag(psi(1,vParaDirichletAlpha));
    
end % end of function


function r = drchrnd(vAlpha,n)
% take a sample from a dirichlet distribution
    p = length(vAlpha);
    r = gamrnd(repmat(vAlpha,n,1),1,n,p);
    r = r ./ repmat(sum(r,2),1,p);

end